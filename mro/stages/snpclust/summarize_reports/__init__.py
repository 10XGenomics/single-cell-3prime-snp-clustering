#!/usr/bin/env python
#
# Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
#

import cellranger.report as cr_report

__MRO__ = """
stage SUMMARIZE_REPORTS(
    in  json cluster_summary,
    in  json evaluate_summary,
    out json summary,
    src py   "stages/snpclust/summarize_reports_pd",
)
"""

def main(args, outs):
    summary_files = [
        args.cluster_summary,
        args.evaluate_summary,
    ]

    cr_report.merge_jsons(summary_files, outs.summary)
