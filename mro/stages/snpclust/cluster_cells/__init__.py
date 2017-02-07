#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import gc
import h5py
import itertools
import json
import numpy as np
import resource
import scipy.stats as sp_stats
import sys
from sklearn.utils.extmath import logsumexp
import tenkit.safe_json as tk_safe_json
import tenkit.stats as tk_stats
import cellranger.matrix as cr_matrix
import snpclust.constants as snp_constants

__MRO__ = '''
stage CLUSTER_CELLS(
    in  path  reference_path,
    in  h5    raw_allele_bc_matrices_h5,
    in  h5    likelihood_allele_bc_matrices_h5,
    in  int   iterations                        "Number of sampling iterations",
    in  int   seed                              "Random seed",
    in  int   burn                              "Number of burn-in iterations",
    in  int   skip                              "Number of iterations to skip for thinning samples",
    in  int   min_depth                         "Minimum depth to use a locus",
    in  float min_posterior                     "Per-cell probability threshold",
    out json  summary,
    src py    "stages/snpclust/cluster_cells_pd",
) split using (
    in  int   seed,
    in  int   components,
    in  int   chain,
)
'''

# Number of sampling chains to run in parallel
NUM_CHAINS = 20

# Clamp posterior probabilities to this value during sampling
MIN_POSTERIOR_PROB = 1e-4

# Maximum number of iterations to run for relabeling algorithm
MAX_RELABEL_ITERATIONS = 20

def normalize_log_probs(log_probs, log=True):
    probs = np.power(10.0, log_probs)
    probs /= probs.sum()
    return np.log10(probs) if log else probs

def fast_choice_among_3(p):
    assert p.shape[1] == 3
    # Sample a k=3 categorical variable given a matrix of probabilities
    # where each row is a vector of 3 probabilities to sample from.
    # 1000x faster than numpy.random.choice
    p0 = p[:,0]
    p1 = p0 + p[:,1]

    u = np.random.uniform(size=p.shape[0])
    bit0 = np.logical_and(u >= p0, u <= p1)
    bit1 = u >= p1

    return bit0.astype(int) + 2*bit1.astype(int)

def split(args):
    chunks = []
    for chain in xrange(NUM_CHAINS):
        for components in [1, 2]:
            chunk = {
                'seed': int(args.seed) + chain,
                'components': components,
                'model': get_model_name(components),
                'chain': chain,
            }
            chunks.append(chunk)

    return {
        'chunks': chunks,
        'join': {
            '__mem_gb': 64, # > K*K!*n*chains*8/1e9 * 3
        },
    }

""" Compute log likelihood of component assignments z conditional on genomes g.
    Store the result in loglk_z_ik
"""
def compute_z_conditional_log_likelihood(K, g_jk, loglk_z_ik, loglk_rr, loglk_ra, loglk_aa):
    for k in xrange(K):
        gk_is_rr = (g_jk[:,k] == 0).astype(float)
        gk_is_ra = (g_jk[:,k] == 1).astype(float)
        gk_is_aa = (g_jk[:,k] == 2).astype(float)
        loglk_z_ik[k,:] = gk_is_rr * loglk_rr + gk_is_ra * loglk_ra + gk_is_aa * loglk_aa

def sample_mixture_prior(K, n, alpha, pi_min, pi_max):
    assert K == 1 or K == 2
    if K == 1:
        pi_k = np.ones(1)
    elif K == 2:
        # for K>2 need dirichlet sample
        pi_0 = sp_stats.beta.rvs(alpha, alpha, scale=(pi_max-pi_min), loc=pi_min)
        pi_k = np.array([pi_0, 1-pi_0])
    return pi_k

def compute_z_conditional_posterior(pi_k, loglk_z_ik):
    logpi_k = np.log(pi_k)
    logpost_z_ik = loglk_z_ik + logpi_k[:,np.newaxis]
    post_z_ik = np.exp(logpost_z_ik - logsumexp(logpost_z_ik, axis=0)[np.newaxis, :])

    # Cap the posterior probability lower bound
    post_z_ik[post_z_ik == 0] = MIN_POSTERIOR_PROB

    # Renormalize
    post_z_ik = post_z_ik / post_z_ik.sum(axis=0)
    return post_z_ik

def compute_complete_log_likelihood(loglk_z_ik, z_i):
    return loglk_z_ik[z_i, np.arange(loglk_z_ik.shape[1])].sum()

def get_log_likelihood_matrices(likelihood_matrices, loci):
    loglk_rr = likelihood_matrices.matrices[snp_constants.HOM_REF_ALLELE].m[loci, :].tocsr()
    loglk_ra = likelihood_matrices.matrices[snp_constants.HET_ALLELE].m[loci, :].tocsr()
    loglk_aa = likelihood_matrices.matrices[snp_constants.HOM_ALT_ALLELE].m[loci, :].tocsr()
    return (loglk_rr, loglk_ra, loglk_aa)

def get_usable_loci(allele_matrices, min_depth):
    r = allele_matrices.matrices[snp_constants.SNP_REF_BASE_TYPE].m
    a = allele_matrices.matrices[snp_constants.SNP_ALT_BASE_TYPE].m

    # Find loci with high enough depth in this sample
    return np.flatnonzero(np.array((r+a).sum(axis=1) >= min_depth))

def sample_model_posterior(likelihood_matrices, allele_matrices, num_components, iterations, burn, thin, min_depth):
    assert num_components == 1 or num_components == 2

    use_loci = get_usable_loci(allele_matrices, min_depth)
    loglk_rr, loglk_ra, loglk_aa = get_log_likelihood_matrices(likelihood_matrices, use_loci)

    # Number of loci
    m = loglk_rr.shape[0]
    # Number of cells
    n = loglk_rr.shape[1]
    # Number of components
    K = num_components

    print "m,n,K = %d,%d,%d" % (m,n,K)
    sys.stdout.flush()

    # Hyperprior on mixture probs pi_k
    alpha = 0.5
    pi_min = 0
    pi_max = 1
    pi_samples = np.zeros((iterations, K))

    # Log likelihood, posterior, and state of cell-genome assignments
    pi_k = np.zeros(K)
    loglk_z_ik = np.zeros((K, n))
    z_i = np.zeros(n, dtype=int)
    z_indicator_samples = np.zeros((iterations, K, n))
    post_z_ik_samples = np.zeros((iterations, K, n))

    # Genotype prior at locus j
    prior_g_j = np.array([0.25, 0.5, 0.25])
    logprior_g_j = np.log(prior_g_j)

    # Log likelihood, posterior, and state of genome genotypes
    loglk_g_jk = np.zeros((3, m, K))
    logpost_g_jk = np.zeros((3, m, K))
    post_g_jk = np.zeros((3, m, K))
    g_jk = np.zeros((m, K), dtype=int)
    g_indicator_samples = np.zeros((iterations, m, 3, K))

    complete_log_likelihood = np.zeros(iterations)

    # Sample initial assignments from prior
    pi_k = sample_mixture_prior(K, n, alpha, pi_min, pi_max)
    z_i = np.random.choice(K, n, p=pi_k)

    num_kept_samples = 0

    for iter_i in xrange(iterations):
        if iter_i % 1000 == 0:
            print "iter %d" % (iter_i)
            sys.stdout.flush()

        # Compute posterior distribution of genomes 1..K conditional on cell-genome assignments z
        for k in xrange(K):
            z_is_k = (z_i == k).astype(float)[:, np.newaxis]
            loglk_g_jk[0,:,k] = (loglk_rr * z_is_k).flatten()
            loglk_g_jk[1,:,k] = (loglk_ra * z_is_k).flatten()
            loglk_g_jk[2,:,k] = (loglk_aa * z_is_k).flatten()

        for k in xrange(K):
            logpost_g_jk[:,:,k] = loglk_g_jk[:,:,k] + logprior_g_j[:, np.newaxis]
            post_g_jk[:,:,k] = np.exp(logpost_g_jk[:,:,k] - logsumexp(logpost_g_jk[:,:,k], axis=0))

        # Sample genomes from conditional posterior
        for k in xrange(K):
            g_jk[:,k] = fast_choice_among_3(post_g_jk[:,:,k].transpose())

        # Compute posterior distribution of cell-genome assignments z conditional on genomes g
        compute_z_conditional_log_likelihood(K, g_jk, loglk_z_ik, loglk_rr, loglk_ra, loglk_aa)

        # Sample mixture probabilities from prior
        pi_k = sample_mixture_prior(K, n, alpha, pi_min, pi_max)
        post_z_ik = compute_z_conditional_posterior(pi_k, loglk_z_ik)

        # Sample cell-genome assignments from conditional posterior
        if K == 2:
            z_i = (np.random.uniform(size=n) >= post_z_ik[0,:]).astype(int)

        # Record selected samples
        if iter_i >= burn and (thin == 0 or (iter_i % thin) == 0):
            for k in xrange(K):
                z_indicator_samples[num_kept_samples, k, np.flatnonzero(z_i == k)] = 1
                g_indicator_samples[num_kept_samples, np.arange(m), g_jk[:,k], k] = 1
            pi_samples[num_kept_samples, :] = pi_k
            post_z_ik_samples[num_kept_samples, :, :] = post_z_ik
            num_kept_samples += 1

        # Record log likelihood for all iterations
        complete_log_likelihood[iter_i] = compute_complete_log_likelihood(loglk_z_ik, z_i)

    return (
        num_kept_samples,
        g_indicator_samples[0:num_kept_samples, :, :, :],
        pi_samples[0:num_kept_samples, :],
        z_indicator_samples[0:num_kept_samples, :],
        post_z_ik_samples[0:num_kept_samples, :, :],  # keep conditional posterior of z for relabelling algo
        complete_log_likelihood, # keep all for log lk
    )

def main(args, outs):
    np.random.seed(args.seed)

    print "Loading matrices"
    sys.stdout.flush()
    allele_matrices = cr_matrix.GeneBCMatrices.load_h5(args.raw_allele_bc_matrices_h5)
    likelihood_matrices = cr_matrix.GeneBCMatrices.load_h5(args.likelihood_allele_bc_matrices_h5)

    f = h5py.File(outs.summary, 'w')
    print "Sampling"
    sys.stdout.flush()
    num_kept_samples, \
        g_indicator_samples, pi_samples, z_indicator_samples, post_z_ik_samples, \
        complete_log_likelihood = sample_model_posterior(likelihood_matrices,
                                                         allele_matrices,
                                                         int(args.components),
                                                         int(args.iterations),
                                                         int(args.burn),
                                                         int(args.skip),
                                                         int(args.min_depth))

    f.create_dataset('chain', data=int(args.chain))
    f.create_dataset('samples', data=num_kept_samples)
    f.create_dataset('g_indicator_samples', data=g_indicator_samples)
    f.create_dataset('pi_samples', data=pi_samples)
    f.create_dataset('z_indicator_samples', data=z_indicator_samples)
    f.create_dataset('post_z_ik_samples', data=post_z_ik_samples)
    f.create_dataset('complete_log_likelihood', data=complete_log_likelihood)
    f.close()

def evaluate_snp_cluster_calls(called, prob, ground_truth, threshold):
    actual = np.array(ground_truth['cluster']).astype(int)
    minor_called_class = 1 - sp_stats.mode(called).mode[0]
    minor_actual_class = 1 - sp_stats.mode(actual).mode[0]
    called = (called == minor_called_class) & (prob >= threshold)
    actual = actual == minor_actual_class

    tp = sum(called & actual)
    tn = sum(np.logical_not(called) & np.logical_not(actual))
    fp = sum(called & np.logical_not(actual))
    fn = sum(np.logical_not(called) & actual)

    return {
        'tp': tp,
        'tn': tn,
        'fp': fp,
        'fn': fn,
        'sensitivity': tk_stats.robust_divide(tp, tp+fn),
        'ppv': tk_stats.robust_divide(tp, tp+fp),
    }

def get_model_name(n_components):
    return 'model%d' % n_components

def summarize_posterior(args, model_name, chunk_defs, chunk_outs, loglk_rr, loglk_ra, loglk_aa, permutations):
    # Summarize the sampled parameter values
    chunks = []
    for chunk_idx, (chunk_def, chunk_out) in enumerate(itertools.izip(chunk_defs, chunk_outs)):
        chunks.append({
            'data': h5py.File(chunk_out.summary, 'r'),
        })

    K = chunks[0]['data']['z_indicator_samples'].shape[1] # components
    n = chunks[0]['data']['z_indicator_samples'].shape[2] # cells
    m = chunks[0]['data']['g_indicator_samples'].shape[1] # loci
    samples_per_chunk = chunks[0]['data']['samples'][()]
    num_samples = samples_per_chunk * len(chunks)
    assert K <= 2
    num_perms = np.math.factorial(K) # 0 = identity permutation, 1 = swapped from original

    z_indicator_sum = np.zeros(chunks[0]['data']['z_indicator_samples'].shape[1:])
    g_indicator_sum = np.zeros(chunks[0]['data']['g_indicator_samples'].shape[1:])

    if permutations is None:
        permutations = np.zeros(num_samples, dtype=int)

    # While summarizing, relabel each chunk according to the given label permutations
    sample_start = 0
    z_indicator_samples_perm = np.zeros((num_perms, samples_per_chunk, K, n), dtype=int)
    g_indicator_samples_perm = np.zeros((num_perms, samples_per_chunk, m, 3, K), dtype=int)
    for chunk in chunks:
        chunk_permutations = permutations[sample_start:(sample_start + samples_per_chunk)]

        z_indicator_samples_perm[0,:,:,:] = chunk['data']['z_indicator_samples'][:]
        if K == 2:
            z_indicator_samples_perm[1,:,:,:] = z_indicator_samples_perm[0,:,::-1,:]
        z_indicator_sum += z_indicator_samples_perm[chunk_permutations,
                                                    np.arange(z_indicator_samples_perm.shape[1]), :, :].sum(axis=0)

        g_indicator_samples_perm[0,:,:,:,:] = chunk['data']['g_indicator_samples'][:]
        if K == 2:
            g_indicator_samples_perm[1,:,:,:,:] = g_indicator_samples_perm[0,:,:,:,::-1]
        g_indicator_sum += g_indicator_samples_perm[chunk_permutations,
                                                    np.arange(g_indicator_samples_perm.shape[1]), :, :, :].sum(axis=0)
        sample_start += samples_per_chunk
        chunk['data'].close()

    # Summarize z
    post_z_ik = z_indicator_sum / z_indicator_sum.sum(axis=0)
    call_z = np.argmax(post_z_ik, axis=0)  # MAP estimate of z
    prob_z = np.max(post_z_ik, axis=0)

    # Summarize g
    post_g_jkg = g_indicator_sum / g_indicator_sum.sum(axis=1)[:,np.newaxis,:]
    call_g = np.argmax(post_g_jkg, axis=1) # MAP estimates of g_k
    prob_g = np.max(post_g_jkg, axis=1)

    # Threshold call on probability
    thresholded_call = [None]*len(call_z)
    for i, p in enumerate(prob_z):
        if prob_z[i] >= args.min_posterior:
            thresholded_call[i] = call_z[i]

    # Summarize pi
    pi_k = np.zeros(K)
    for k in xrange(K):
        pi_k[k] = np.mean(call_z == k)

    d = {}
    d.update({
        'probability': prob_z.astype(float).tolist(),
        'call': call_z.astype(int).tolist(),
        'thresholded_call': thresholded_call,
    })
    for k in xrange(K):
        d.update({('g%d_probability' % k): prob_g[:,k].astype(float).tolist(),
                  ('g%d_call' % k): call_g[:,k].astype(int).tolist(),
              })
    d['pi'] = pi_k.astype(float).tolist()

    # Mix fraction (fraction of minor population)
    mix_fraction_all_calls = tk_stats.robust_divide(np.sum(call_z == 0), len(call_z))
    mix_fraction_all_calls = min(mix_fraction_all_calls, 1-mix_fraction_all_calls)
    mix_fraction = tk_stats.robust_divide(np.sum(np.logical_and(call_z == 0, prob_z >= args.min_posterior)), np.sum(prob_z >= args.min_posterior))
    mix_fraction = min(mix_fraction, 1-mix_fraction)
    d['mix_fraction_all_calls'] = mix_fraction_all_calls
    d['mix_fraction'] = mix_fraction

    d['relabel_permutation'] = permutations.astype(int).tolist()

    # Compute message length
    loglk_z_ik = np.zeros((K, n))
    compute_z_conditional_log_likelihood(K, call_g, loglk_z_ik, loglk_rr, loglk_ra, loglk_aa)
    d['complete_log_likelihood'] = compute_complete_log_likelihood(loglk_z_ik, call_z)
    d['message_length'] = compute_message_length(d['complete_log_likelihood'], n, pi_k, m)

    return {("%s_%s" % (model_name, k)):v for k,v in d.iteritems()}

def compute_message_length(log_likelihood, n_obs, mixing_probs, component_dimension):
    # Message length for a categorical mixture model
    # Formula given in http://arxiv.org/pdf/1409.7419.pdf
    M = float(component_dimension)
    n = float(n_obs)
    alpha = mixing_probs
    k_nz = np.sum(mixing_probs > 0)

    alpha_sum = 0
    for k in np.flatnonzero(alpha > 0):
        alpha_sum += np.log(n*alpha[k]/12)

    return (M/2.0)*alpha_sum + (k_nz/2.0)*np.log(n/12.0) + (k_nz*(M+1.0))/2.0 - log_likelihood

def relabel_samples(chunk_defs, chunk_outs):
    # Relabel mixture model samples by minimizing the K-L divergence from the posterior mean
    # Algorithm described in http://stephenslab.uchicago.edu/assets/papers/Stephens2000b.pdf

    # Get sampled z conditional posteriors
    chunks = []
    for chunk_idx, (chunk_def, chunk_out) in enumerate(itertools.izip(chunk_defs, chunk_outs)):
        if chunk_def.components == 2:
            chunks.append({
                'data': h5py.File(chunk_out.summary, 'r'),
                'chunk_idx': chunk_idx,
                'chunk_def': chunk_def,
                'chunk_outs': chunk_out,
            })

    _, K, n = chunks[0]['data']['z_indicator_samples'].shape
    assert K == 2
    samples_per_chunk = chunks[0]['data']['samples'][()]
    num_samples = sum(chunk['data']['samples'][()] for chunk in chunks)
    num_perms = np.math.factorial(K) # 0 = identity permutation, 1 = swapped from original

    # Load data from all chunks
    gc.collect()
    post_z_ik_samples_perm = np.zeros((num_perms, num_samples, K, n))
    print "Loading all chunk data"
    print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    sys.stdout.flush()

    sample_start = 0
    for i, chunk in enumerate(chunks):
        chunk_post_z_ik_samples = chunk['data']['post_z_ik_samples'][:]
        post_z_ik_samples_perm[0,sample_start:(sample_start+samples_per_chunk),:,:] = chunk_post_z_ik_samples
        post_z_ik_samples_perm[1,sample_start:(sample_start+samples_per_chunk),:,:] = chunk_post_z_ik_samples[:,::-1,:]
        del chunk_post_z_ik_samples
        gc.collect()
        sample_start += samples_per_chunk
        print "%d\t%d" % (i, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        sys.stdout.flush()

    print "done"
    print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    sys.stdout.flush()

    for chunk in chunks:
        chunk['data'].close()

    # Initialize w/ the identity permutation
    permutations = np.zeros(num_samples, dtype=int)
    prev_permutations = np.zeros(num_samples, dtype=int)
    converged = False
    iterations = 0

    for _ in xrange(MAX_RELABEL_ITERATIONS):
        # Step 1: Compute q_hat, the posterior mean classification probabilities, given the current permutations
        q_hat_ik = np.mean(post_z_ik_samples_perm[permutations,np.arange(num_samples),:,:], axis=0)
        sys.stdout.write('.')
        sys.stdout.flush()

        # Step 2: Find the permutations that minimize the K-L divergence of the sampled parameters from their posterior means
        kl_div = (post_z_ik_samples_perm * np.log(post_z_ik_samples_perm / q_hat_ik)).sum(axis=(2,3))
        permutations = np.argmin(kl_div, axis=0)
        sys.stdout.write(',')
        sys.stdout.flush()

        iterations += 1
        permutations_changed = np.sum(permutations != prev_permutations)
        print 'relabel iter %d; changed %d; sum(qhat)=%0.4f, kl_div=%0.4f' % (iterations, permutations_changed, q_hat_ik.sum(), kl_div[permutations,np.arange(num_samples)].sum())

        del kl_div
        del q_hat_ik
        gc.collect()
        print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        sys.stdout.flush()

        if np.all(permutations == prev_permutations):
            print 'Converged after %d' % iterations
            converged = True
            break
        np.copyto(prev_permutations, permutations)
    if not converged:
        print 'Warning: relabelling did not converge after %d iterations' % iterations

    return permutations, converged, iterations

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    allele_matrices = cr_matrix.GeneBCMatrices.load_h5(args.raw_allele_bc_matrices_h5)
    likelihood_matrices = cr_matrix.GeneBCMatrices.load_h5(args.likelihood_allele_bc_matrices_h5)

    d = {}
    models = [get_model_name(K) for K in [1,2]]

    for model in models:
        print "Summarizing %s" % (model)
        sys.stdout.flush()
        use_loci = get_usable_loci(allele_matrices, int(args.min_depth))
        loglk_rr, loglk_ra, loglk_aa = get_log_likelihood_matrices(likelihood_matrices, use_loci)
        chunks = [(chunk_def, chunk_out) for chunk_def, chunk_out in itertools.izip(chunk_defs, chunk_outs) if chunk_def.model == model]

        # Concatenate the chains and relabel the full set of samples
        if chunks[0][0].components > 1:
            permutations, converged, iterations = relabel_samples([chunk[0] for chunk in chunks],
                                                                  [chunk[1] for chunk in chunks])
            d['%s_relabel_converged' % model] = converged
            d['%s_relabel_iterations' %  model] = iterations
        else:
            permutations = None

        # Generate the summary
        d.update(summarize_posterior(args, model,
                                     [chunk[0] for chunk in chunks],
                                     [chunk[1] for chunk in chunks],
                                     loglk_rr, loglk_ra, loglk_aa,
                                     permutations
                                 ))

    # Select model based on minimum message length
    message_lengths = np.array([d['%s_message_length' % model_name] for model_name in models], dtype=float)
    d['selected_model'] = models[np.argmin(message_lengths)]

    with open(outs.summary, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(d), f, indent=4, sort_keys=True)
