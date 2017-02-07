#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#
import json
import numpy as np
import scipy.stats as sp_stats
import tenkit.stats as tk_stats

__MRO__ = '''
stage EVALUATE_SNP_CLUSTERS(
    in  json cluster_summary,
    in  json ground_truth,
    out json summary,
    src py   "stages/snpclust/evaluate_snp_clusters_pd",
)
'''

def evaluate_snp_cluster_calls(cluster_assignment, thresholded_calls, actual):
    """ Args:
        - cluster_assignment: list(int)
        - thresholded_calls: list(int), None if no call
        - actual: list(int) """
    cluster_assignment = np.array(cluster_assignment, dtype=int)
    actual = np.array(actual, dtype=int)

    minor_called_class = 1 - sp_stats.mode(cluster_assignment).mode[0]
    minor_actual_class = 1 - sp_stats.mode(actual).mode[0]

    was_called = np.array([x is not None for x in thresholded_calls])

    called_pos = (cluster_assignment == minor_called_class)[was_called]
    actual_pos = (actual == minor_actual_class)[was_called]

    nc = sum(np.logical_not(was_called))
    tp = sum(called_pos & actual_pos)
    tn = sum(np.logical_not(called_pos) & np.logical_not(actual_pos))
    fp = sum(called_pos & np.logical_not(actual_pos))
    fn = sum(np.logical_not(called_pos) & actual_pos)

    return {
        'tp': tp,
        'tn': tn,
        'fp': fp,
        'fn': fn,
        'sensitivity': tk_stats.robust_divide(tp, tp+fn),
        'ppv': tk_stats.robust_divide(tp, tp+fp),
        'no_call_rate': tk_stats.robust_divide(nc, len(actual)),
    }

def main(args, outs):
    if args.ground_truth is None:
        # Write empty json
        with open(outs.summary, 'w') as f:
            f.write('{}')
        return

    with open(args.cluster_summary) as f:
        cluster_summary = json.load(f)

    with open(args.ground_truth) as f:
        ground_truth = json.load(f)['cluster']

    summary = {}

    for model in ['model1', 'model2']:
        model_summary = evaluate_snp_cluster_calls(cluster_summary['%s_call' % model],
                                                   cluster_summary['%s_thresholded_call' % model],
                                                   ground_truth)
        summary.update({("%s_%s" % (model, key)):val for key,val in model_summary.iteritems()})

    with open(outs.summary, 'w') as f:
        json.dump(summary, f, indent=4, sort_keys=True)
