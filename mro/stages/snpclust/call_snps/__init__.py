#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import martian
import subprocess
import tenkit.bam as tk_bam
import tenkit.bio_io as tk_io
import tenkit.constants as tk_constants
import cellranger.utils as cr_utils

__MRO__ = '''
stage CALL_SNPS(
    in  path   reference_path,
    in  bam    input,
    in  int    n_donors,
    out vcf[]  output,
    src py     "stages/snpclust/call_snps_pd",
) split using (
    in  string locus,
)
'''

def split(args):
    in_bam = tk_bam.create_bam_infile(args.input)
    loci = tk_bam.generate_tiling_windows(in_bam, tk_constants.PARALLEL_LOCUS_SIZE)
    chunks = [{'locus': locus} for locus in loci]
    return {'chunks': chunks}

def main(args, outs):
    genome_fasta_path = cr_utils.get_reference_genome_fasta(args.reference_path)

    chrom, start, stop = tk_io.get_locus_info(args.locus)
    bed_path = martian.make_path('region.bed')
    with open(bed_path, 'w') as f:
        f.write(chrom+"\t"+str(start)+"\t"+str(stop)+"\n")

    freebayes_args = ['freebayes', '-f', genome_fasta_path, '-b', args.input, '-0', '-t', bed_path]

    with open(outs.output, 'w') as f:
        subprocess.check_call(freebayes_args, stdout=f)

def join(args, outs, chunk_defs, chunk_outs):
    outs.output = [chunk.output for chunk in chunk_outs]
