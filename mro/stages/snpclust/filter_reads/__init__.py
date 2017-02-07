#!/usr/bin/env python
#
# Copyright (c) 2016 10X Genomics, Inc. All rights reserved.
#
import itertools
import tenkit.bam as tk_bam
import cellranger.constants as cr_constants
import cellranger.utils as cr_utils

__MRO__ = '''
stage FILTER_READS(
    in  bam    input,
    in  tsv    cell_barcodes,
    in  map    align,
    out bam    output,
    src py     "stages/snpclust/filter_reads_pd",
) split using (
    in  string chunk_start,
    in  string chunk_end,
)
'''

def split(args):
    with tk_bam.create_bam_infile(args.input) as in_bam:
        chunks = tk_bam.chunk_bam_records(in_bam, chunk_bound_key=cr_utils.pos_sort_key,
                                          chunk_size_gb=cr_constants.BAM_CHUNK_SIZE_GB,
                                          max_chunks=cr_constants.MAX_BAM_CHUNKS)
    return {'chunks': chunks}

def main(args, outs):
    outs.coerce_strings()

    in_bam = tk_bam.create_bam_infile(args.input)
    in_bam_chunk = tk_bam.read_bam_chunk(in_bam, (args.chunk_start, args.chunk_end))
    out_bam, _ = tk_bam.create_bam_outfile(outs.output, None, None, template=in_bam)
    cell_bcs = set(cr_utils.load_barcode_tsv(args.cell_barcodes))

    for (tid, pos), reads_iter in itertools.groupby(in_bam_chunk, key=cr_utils.pos_sort_key):
        dupe_keys = set()
        for read in reads_iter:
            if cr_utils.get_read_barcode(read) not in cell_bcs:
                continue

            if cr_utils.is_read_dupe_candidate(read, cr_utils.get_high_conf_mapq(args.align)):
                dupe_key = (cr_utils.si_pcr_dupe_func(read), cr_utils.get_read_umi(read))
                if dupe_key in dupe_keys:
                    continue

                dupe_keys.add(dupe_key)
                out_bam.write(read)

def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    input_bams = [str(chunk.output) for chunk in chunk_outs]
    tk_bam.concatenate(outs.output, input_bams)
    tk_bam.index(outs.output)
