#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import collections
import json
import martian
import numpy as np
import scipy.misc as sp_misc
import scipy.stats as sp_stats
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.constants as tk_constants
import cellranger.constants as cr_constants
import cellranger.matrix as cr_matrix
import cellranger.utils as cr_utils
import snpclust.constants as snp_constants

__MRO__ = '''
stage COUNT_ALLELES(
    in  path  reference_path,
    in  bam   reads,
    in  vcf[] variants,
    in  tsv   cell_barcodes,
    in  int   min_snp_call_qual,
    in  int   min_bcs_per_snp,
    in  int   min_snp_obs,
    in  int   min_snp_base_qual,
    in  float base_error_rate,
    out vcf   filtered_variants,
    out h5    raw_allele_bc_matrices_h5,
    out path  raw_allele_bc_matrices_mex,
    out h5    likelihood_allele_bc_matrices_h5,
    out path  likelihood_allele_bc_matrices_mex,
    src py    "stages/snpclust/count_alleles_pd",
) split using (
    in  vcf   chunk_variants,
    in  json  snps,
)
'''

def format_record(record):
    return ','.join((record.CHROM, str(record.POS), str(record.REF), str(record.ALT[0])))

def vcf_record_iter(in_filename, min_snp_qual):
    in_vcf = tk_io.VariantFileReader(in_filename)
    for record in in_vcf.record_getter(restrict_type='snp'):
        # Only support 1 ALT
        if len(record.ALT) > 1:
            continue
        assert len(record.ALT) == 1

        # Filter SNP based on call quality
        if record.QUAL < min_snp_qual:
            continue

        yield record

def split(args):
    out_json_file = martian.make_path('snps.json')
    min_snp_call_qual = args.min_snp_call_qual if args.min_snp_call_qual is not None \
                        else snp_constants.DEFAULT_MIN_SNP_CALL_QUAL
    save_snps(out_json_file, args.variants, min_snp_call_qual)

    chunks = [{'chunk_variants': chunk_variants, 'snps': out_json_file} for chunk_variants in args.variants]
    return {'chunks': chunks}

def save_snps(out_filename, in_filenames, min_snp_call_qual):
    snps = []
    for in_filename in in_filenames:
        for record in vcf_record_iter(in_filename, min_snp_call_qual):
            snps.append(format_record(record))

    with open(out_filename, 'w') as f:
        json.dump(snps, f)

def load_snps(filename):
    # HACK: Save SNPs as Gene tuples so we can reuse code in GeneBCMatrices
    with open(filename, 'r') as f:
        return [cr_constants.Gene(str(snp), '', None, None, None) for snp in json.load(f)]

def get_read_qpos(read):
    # Support pysam 0.7.8 and 0.9.0
    if hasattr(read, 'qpos'):
        return read.qpos
    elif hasattr(read, 'query_position'):
        return read.query_position

    raise Exception("Pysam pileup read has neither qpos nor query_position attribute")

def main(args, outs):
    in_bam = tk_bam.create_bam_infile(args.reads)

    out_vcf = tk_io.VariantFileWriter(
        open(outs.filtered_variants, 'w'), template_file=open(args.chunk_variants))

    snps = load_snps(args.snps)
    bcs = cr_utils.load_barcode_tsv(args.cell_barcodes)

    raw_matrix_types = snp_constants.SNP_BASE_TYPES
    raw_matrix_snps = [snps for _ in snp_constants.SNP_BASE_TYPES]
    raw_allele_bc_matrices = cr_matrix.GeneBCMatrices(raw_matrix_types, raw_matrix_snps, bcs)

    likelihood_matrix_types = snp_constants.ALLELES
    likelihood_matrix_snps = [snps for _ in snp_constants.ALLELES]
    likelihood_allele_bc_matrices = cr_matrix.GeneBCMatrices(likelihood_matrix_types, likelihood_matrix_snps, bcs, dtype=np.float64)

    # Configurable SNP filter parameters
    min_snp_call_qual = args.min_snp_call_qual if args.min_snp_call_qual is not None else snp_constants.DEFAULT_MIN_SNP_CALL_QUAL
    min_bcs_per_snp = args.min_bcs_per_snp if args.min_bcs_per_snp is not None else snp_constants.DEFAULT_MIN_BCS_PER_SNP
    min_snp_obs = args.min_snp_obs if args.min_snp_obs is not None else snp_constants.DEFAULT_MIN_SNP_OBS
    base_error_rate = args.base_error_rate if args.base_error_rate is not None else snp_constants.DEFAULT_BASE_ERROR_RATE
    min_snp_base_qual = args.min_snp_base_qual if args.min_snp_base_qual is not None else snp_constants.DEFAULT_MIN_SNP_BASE_QUAL

    for record in vcf_record_iter(args.chunk_variants, min_snp_call_qual):
        ref_base = str(record.REF)
        alt_base = str(record.ALT[0])

        pos = record.POS - 1
        snps = collections.defaultdict(lambda: np.zeros((2, 2)))
        for col in in_bam.pileup(record.CHROM, pos, pos+1):
            if col.pos != pos:
                continue

            for read in col.pileups:
                bc = cr_utils.get_read_barcode(read.alignment)
                umi = cr_utils.get_read_umi(read.alignment)
                assert bc in set(bcs) and umi is not None

                # Overlaps an exon junction
                qpos = get_read_qpos(read)
                if qpos is None:
                    continue

                base = str(read.alignment.query[qpos - read.alignment.qstart])
                base_qual = ord(read.alignment.qual[qpos - read.alignment.qstart]) - tk_constants.ILLUMINA_QUAL_OFFSET

                if base == ref_base:
                    base_index = 0
                elif base == alt_base:
                    base_index = 1
                else:
                    continue

                dupe_key = (bc, umi)
                snps[dupe_key][base_index, 0] += 1
                snps[dupe_key][base_index, 1] = max(base_qual, snps[dupe_key][base_index, 1])

        bcs_bases = collections.defaultdict(collections.Counter)
        for (bc, umi), bases in snps.iteritems():
            base_index = np.argmax(bases[:, 0])
            base = ref_base if base_index == 0 else alt_base
            base_qual = bases[base_index, 1]
            if base_qual < min_snp_base_qual:
                continue
            bcs_bases[bc][base] += 1

        # Filter if not enough unique barcodes
        if len(bcs_bases) < min_bcs_per_snp:
            continue

        # Filter if not enough observed bases
        snp_obs = 0
        for b in bcs_bases.itervalues():
            snp_obs += sum([count for count in b.itervalues()])
        if snp_obs < min_snp_obs:
            continue

        for bc, bases in bcs_bases.iteritems():
            ref_obs = bases[ref_base]
            alt_obs = bases[alt_base]
            total_obs = ref_obs + alt_obs
            obs = np.array([
                ref_obs,
                alt_obs,
            ])

            log_p_hom_ref = sp_stats.binom.logpmf(ref_obs, total_obs, 1 - base_error_rate)
            log_p_hom_alt = sp_stats.binom.logpmf(alt_obs, total_obs, 1 - base_error_rate)
            log_p_het = sp_stats.binom.logpmf(ref_obs, total_obs, 0.5)

            log_p = np.array([
                log_p_hom_ref,
                log_p_het,
                log_p_hom_alt,
            ])
            log_p -= sp_misc.logsumexp(log_p)

            matrix = raw_allele_bc_matrices.matrices.values()[0]
            snp_index = matrix.gene_id_to_int(format_record(record))
            bc_index = matrix.bc_to_int(bc)

            for i, base_type in enumerate(snp_constants.SNP_BASE_TYPES):
                raw_allele_bc_matrices.get_matrix(base_type).m[snp_index, bc_index] = obs[i]

            for i, allele in enumerate(snp_constants.ALLELES):
                likelihood_allele_bc_matrices.get_matrix(allele).m[snp_index, bc_index] = log_p[i]

        out_vcf.write_record(record)

    raw_allele_bc_matrices.save_h5(outs.raw_allele_bc_matrices_h5)
    likelihood_allele_bc_matrices.save_h5(outs.likelihood_allele_bc_matrices_h5)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()

    input_vcfs = [chunk_out.filtered_variants for chunk_out in chunk_outs]
    tk_io.combine_vcfs(outs.filtered_variants, input_vcfs)

    raw_chunk_h5s = [chunk_out.raw_allele_bc_matrices_h5 for chunk_out in chunk_outs]
    raw_allele_bc_matrices = cr_matrix.merge_matrices(raw_chunk_h5s)

    likelihood_chunk_h5s = [chunk_out.likelihood_allele_bc_matrices_h5 for chunk_out in chunk_outs]
    likelihood_allele_bc_matrices = cr_matrix.merge_matrices(likelihood_chunk_h5s)

    raw_allele_bc_matrices.save_h5(outs.raw_allele_bc_matrices_h5)
    raw_allele_bc_matrices.save_mex(outs.raw_allele_bc_matrices_mex)
    likelihood_allele_bc_matrices.save_h5(outs.likelihood_allele_bc_matrices_h5)
    likelihood_allele_bc_matrices.save_mex(outs.likelihood_allele_bc_matrices_mex)
