#!/usr/bin/env python
#
# Copyright (c) 2015 10X Genomics, Inc. All rights reserved.
#
import collections
import cPickle
import h5py as h5
import itertools
import json
import math
import numpy as np
import operator
import random
import tenkit.fasta as tk_fasta
import tenkit.safe_json as tk_safe_json
import tenkit.seq as tk_seq
import tenkit.stats as tk_stats
import cellranger.constants as cr_constants
import cellranger.stats as cr_stats
import cellranger.utils as cr_utils

# Metrics
class Metric:
    def __init__(self, report_type=cr_constants.DEFAULT_REPORT_TYPE, **kwargs):
        self.report_type = report_type
        self.active = kwargs.get('always_active', False)

    def add(self, elem):
        self.active = True

    def merge(self, metric):
        raise NotImplementedError

    def report(self):
        raise NotImplementedError

# Dictionary metrics i.e. histograms, element counts
class DictionaryMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.d = collections.Counter()

    def add(self, elem, value=1):
        Metric.add(self, elem)
        self.d[elem] += value

    def merge(self, metric):
        for elem in metric.d:
            self.d[elem] += metric.d[elem]

    def report(self):
        return self.d

class PercentDictionaryMetric(DictionaryMetric):
    def report(self):
        total = sum(self.d.values())

        d = {}
        for key, value in self.d.iteritems():
            d[key] = tk_stats.robust_divide(float(value), float(total))
        return d

class HistogramMetric(DictionaryMetric):
    def __init__(self, cutoffs, **kwargs):
        DictionaryMetric.__init__(self, **kwargs)
        self.cutoffs = cutoffs
        for cutoff in self.cutoffs:
            self.d[cutoff] = 0

    def add(self, elem, value=1):
        Metric.add(self, elem)
        for cutoff in reversed(self.cutoffs):
            if elem >= cutoff:
                self.d[cutoff] += value
                break

class NumpyHistogramMetric(Metric):
    """ A histogram with bins of size 1 and a final bin containing elements > max_value """
    def __init__(self, max_value=10, **kwargs):
        Metric.__init__(self, **kwargs)
        self.counts = np.zeros(2+max_value, dtype=np.int_)
        self.max_value = max_value

    def add(self, elem):
        assert elem >= 0
        Metric.add(self, elem)
        if elem > self.max_value:
            self.counts[-1] += 1
        else:
            self.counts[elem] += 1

    def add_many(self, elems):
        self.active = True
        elems = np.copy(elems).astype(np.int_)
        elems[elems > self.max_value] = 1 + self.max_value
        self.counts += np.bincount(elems, minlength=len(self.counts))

    def merge(self, metric):
        assert self.max_value == metric.max_value
        self.counts += metric.counts

    def report(self):
        d = {str(k):int(v) for k, v in itertools.izip(xrange(0, 1 + self.max_value), self.counts)}
        d[">%d" % self.max_value] = int(self.counts[-1])
        return d

class SequenceDistributionMetric(Metric):
    """ Count base frequencies at k positions """
    def __init__(self, k, **kwargs):
        Metric.__init__(self, **kwargs)
        assert k > 0
        self.k = k
        self.counts = np.zeros((len(tk_seq.NUCS), self.k), dtype=np.int_)

    def add(self, kmer, value=1):
        Metric.add(self, kmer)
        self.counts[[tk_seq.NUCS_INVERSE[nuc] for nuc in kmer], np.arange(0, self.k)] += int(value)

    def merge(self, metric):
        self.counts += metric.counts

    def report(self):
        d = collections.defaultdict(lambda: collections.defaultdict(int))
        for col, pos in enumerate(xrange(-self.k/2, self.k/2)):
            for nuc, row in tk_seq.NUCS_INVERSE.iteritems():
                d[str(pos)][nuc] = self.counts[row, col]
        return d

class MedianMetric(DictionaryMetric):
    def report(self):
        return cr_stats.compute_median_from_distribution(self.d)

class IQRMetric(DictionaryMetric):
    def report(self):
        return cr_stats.compute_iqr_from_distribution(self.d)

class TopNMetric(DictionaryMetric):
    def __init__(self, topN, **kwargs):
        DictionaryMetric.__init__(self, **kwargs)
        self.topN = topN

    def report(self):
        items = sorted(self.d.items(), key=operator.itemgetter(1), reverse=True)
        topN = {key: value for key, value in items[:self.topN]}
        return topN

class SubsampledTopNMetric(TopNMetric):
    def __init__(self, topN, sample_rate, **kwargs):
        TopNMetric.__init__(self, topN, **kwargs)
        self.sample_rate = sample_rate

    def add(self, elem, value=1):
        if cr_utils.downsample(self.sample_rate):
            TopNMetric.add(self, elem, value)

class EffectiveDiversityMetric(DictionaryMetric):
    def report(self):
        counts = np.array(self.d.values())
        return cr_stats.effective_diversity(counts)

class SetMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.d = set()

    def add(self, elem):
        Metric.add(self, elem)
        self.d.add(elem)

    def merge(self, metric):
        self.d.update(metric.d)

    def report(self):
        return len(self.d)

# Stats metrics i.e. mean, total count, percent, sum, stddev
class StatsMetric(Metric):
    def __init__(self, **kwargs):
        Metric.__init__(self, **kwargs)
        self.m1 = 0
        self.m0 = 0

    def merge(self, metric):
        self.m1 += metric.m1
        self.m0 += metric.m0

class PercentMetric(StatsMetric):
    def add(self, elem, filter=False):
        Metric.add(self, elem)
        self.m0 += elem
        if filter:
            self.m1 += elem

    def set_value(self, numerator, denominator):
        self.active = True
        self.m0 = denominator
        self.m1 = numerator

    def report(self):
        return float(self.m1) / float(self.m0) if self.m0 > 0 else 0

class RateMetric(StatsMetric):
    def add(self, elem, numerator=False):
        Metric.add(self, elem)
        if numerator:
            self.m1 += elem
        else:
            self.m0 += elem

    def report(self):
        return float(self.m1) / float(self.m0) if self.m0 > 0 else 0

class CountMetric(StatsMetric):
    def add(self, elem):
        Metric.add(self, elem)
        self.m0 += 1

    def set_value(self, value):
        self.active = True
        self.m0 = value

    def report(self):
        return self.m0

class MeanMetric(StatsMetric):
    def add(self, elem, weight=1):
        Metric.add(self, elem)
        self.m1 += elem * weight
        self.m0 += weight

    def add_many(self, elems):
        self.active = True
        self.m1 += elems.sum()
        self.m0 += len(elems)

    def report(self):
        return float(self.m1) / float(self.m0) if self.m0 > 0 else 0.0

class VarianceMetric(StatsMetric):
    """ Chan et. al (1979) parallel variance """
    def __init__(self, **kwargs):
        StatsMetric.__init__(self, **kwargs)
        self.m2 = 0

    def add(self, elem):
        Metric.add(self, elem)
        self.m0 += 1
        delta = float(elem - self.m1)
        self.m1 += delta / float(self.m0)
        self.m2 += delta * (float(elem - self.m1))

    def add_many(self, elems):
        self.active = True
        elems_mean = elems.mean()
        elems_m2 = np.square(elems - elems_mean).sum()
        delta = elems_mean - self.m1
        self.m1 = (self.m0 * self.m1 + len(elems) * elems_mean) / float(self.m0 + len(elems))
        self.m2 += elems_m2 + delta * delta * float(self.m0) * float(len(elems)) / (self.m0 + len(elems))
        self.m0 += len(elems)

    def merge(self, metric):
        delta = metric.m1 - self.m1
        self.mean = (self.m0 * self.m1 + metric.m0 * metric.m1) / float(self.m0 + metric.m0)
        self.m2 += metric.m2 + delta * delta * float(self.m0) * float(metric.m0) / (self.m0 + metric.m0)
        self.m0 += metric.m0

    def report(self):
        if self.m0 < 2:
            return float('nan')
        else:
            return tk_stats.robust_divide(self.m2, self.m0 - 1)

class CVMetric(VarianceMetric):
    """ Coeffcient of variation: sd/mean """
    def report(self):
        if self.m0 < 2:
            return float('nan')
        else:
            sd = math.sqrt(VarianceMetric.report(self))
            return tk_stats.robust_divide(sd, self.m1)

class DispersionMetric(VarianceMetric):
    """ Method of moments estimator for negative binomial 1/r """
    def report(self):
        if self.m0 < 2:
            return float('nan')
        else:
            var = VarianceMetric.report(self)
            return (var - self.m1) / (self.m1 * self.m1)

# For embedding metrics within one another
class EmbeddedMetric(Metric):
    def __init__(self, metric_cls, **kwargs):
        Metric.__init__(self, **kwargs)
        self.d = {}
        self.metric_cls = metric_cls

    def add(self, elem, *args):
        Metric.add(self, elem)
        if elem not in self.d:
            self.d[elem] = self.metric_cls(report_type=self.report_type)
        self.d[elem].add(*args)

    def merge(self, metric):
        for elem in metric.d:
            if elem not in self.d:
                self.d[elem] = metric.d[elem]
            else:
                self.d[elem].merge(metric.d[elem])

    def report(self):
        return {key: metric.report() for key, metric in self.d.iteritems()}

# For metrics binned by a numeric value
class BinnedMetric(EmbeddedMetric):
    def __init__(self, metric_cls, cutoffs, **kwargs):
        EmbeddedMetric.__init__(self, metric_cls, **kwargs)
        self.cutoffs = cutoffs
        for cutoff in cutoffs:
            self.d[cutoff] = self.metric_cls(report_type=self.report_type)

    def add(self, elem, *args):
        Metric.add(self, elem)
        for cutoff in reversed(self.cutoffs):
            if elem >= cutoff:
                self.d[cutoff].add(*args)
                break

METRICS = {
    # Raw fastq metrics
    'total_reads':               (CountMetric, {}),
    'total_reads_per_gem_group': (CountMetric, {'prefixes': ['gem_groups']}),
    'total_read_pairs':          (CountMetric, {}),
    'read_bases_with_q30_frac':  (PercentMetric, {}),
    'read_N_bases_frac':         (PercentMetric, {}),
    'top_read_prefixes':         (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),
    'read2_bases_with_q30_frac': (PercentMetric, {}),
    'read2_N_bases_frac':        (PercentMetric, {}),
    'top_read2_prefixes':        (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),

    # Primer fastq metrics
    'read_with_perfect_primer_or_homopolymer_frac':  (PercentMetric, {}),
    'perfect_homopolymer_frac':                      (PercentMetric, {'prefixes': ['nucs']}),
    'perfect_primer_frac':                           (PercentMetric, {'prefixes': ['primer_names']}),

    # Aligned bam metrics
    'good_bc_frac':            (PercentMetric, {}),
    'barcodes_detected':       (SetMetric, {}),
    'corrected_bc_frac':       (PercentMetric, {}),
    'bc_bases_with_q30_frac':  (PercentMetric, {}),
    'bc_N_bases_frac':         (PercentMetric, {}),

    'sample_index_bases_with_q30_frac':   (PercentMetric, {}),
    'sample_index_N_bases_frac':          (PercentMetric, {}),

    'good_umi_frac':            (PercentMetric, {}),
    'corrected_umi_frac':       (PercentMetric, {}),
    'umi_bases_with_q30_frac':  (PercentMetric, {}),
    'umi_N_bases_frac':         (PercentMetric, {}),
    'failed_rescue_frac':       (PercentMetric, {}),

    # Transcript-based metrics
    'read_distance_to_transcript_three_prime': (HistogramMetric, {'args': [cr_constants.READ_POSITION_CUTOFFS]}),
    'read_distance_to_transcript_five_prime': (HistogramMetric, {'args': [cr_constants.READ_POSITION_CUTOFFS]}),

    # Hierarchical (exclusive) sequence filters
    'umi_filter_frac':            (PercentMetric, {'prefixes': ['umi_properties']}),
    'barcode_filter_frac':        (PercentMetric, {'prefixes': ['barcode_properties']}),

    # Non-hierarchical (non-exclusive) sequence properties
    'umi_property_frac':          (PercentMetric, {'prefixes': ['umi_properties']}),
    'barcode_property_frac':      (PercentMetric, {'prefixes': ['barcode_properties']}),

    'median_insert_size':    (MedianMetric, {'prefixes': ['references']}),
    'insert_size_histogram': (HistogramMetric, {'args': [cr_constants.INSERT_SIZE_CUTOFFS], 'prefixes': ['references']}),
    'iqr_insert_size':       (IQRMetric, {'prefixes': ['references']}),

    'reads_frac': (PercentMetric, {'prefixes': ['references', 'regions', 'read_types']}),
    'unmapped_reads_frac': (PercentMetric, {}),
    'antisense_reads_frac': (PercentMetric, {'prefixes': ['references']}),

    'dupe_reads_frac': (PercentMetric, {'prefixes': ['references', 'dupe_types']}),
    'dupe_reads_rate': (RateMetric,{'prefixes': ['references', 'dupe_types']}),

    'umi_hamming_distance_per_dupe_group_histogram': (DictionaryMetric, {'prefixes': ['references', 'dupe_types']}),
    'reads_per_dupe_group_histogram': (NumpyHistogramMetric, {'args': [cr_constants.PER_DUPE_GROUP_MAX], 'prefixes': ['references', 'dupe_types']}),
    'umis_per_dupe_group_histogram': (NumpyHistogramMetric, {'args': [cr_constants.PER_DUPE_GROUP_MAX], 'prefixes': ['references', 'dupe_types']}),
    'reads_per_molecule_histogram': (NumpyHistogramMetric, {'args': [cr_constants.PER_DUPE_GROUP_MAX], 'prefixes': ['references', 'dupe_types']}),

    'top_raw_umis':       (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),
    'top_processed_umis': (TopNMetric, {'args': [cr_constants.TOP_N], 'prefixes': ['references', 'read_types']}),
    'effective_umi_diversity': (EffectiveDiversityMetric, {'prefixes': ['references', 'read_types']}),

    'top_raw_barcodes':            (SubsampledTopNMetric, {'args': [cr_constants.TOP_N, cr_constants.TOP_RAW_SEQUENCE_SAMPLE_RATE]}),
    'effective_barcode_diversity': (EffectiveDiversityMetric, {'prefixes': ['references', 'read_types']}),

    'barcode_reads': (DictionaryMetric, {'kwargs': {'report_type': 'barcodes', 'always_active': True}, 'prefixes': ['references', 'regions', 'read_types']}),

    ## Subsampled metrics
    'subsampled_duplication_frac': (PercentMetric, {'kwargs': {'always_active': True}, 'prefixes': ['references', 'subsample_types', 'subsample_depths']}),

    # Subsampled matrix summary metrics
    'subsampled_filtered_bcs_median_unique_genes_detected': (CountMetric, {'kwargs': {'always_active': True}, 'prefixes': ['references', 'subsample_types', 'subsample_depths']}),
    'subsampled_filtered_bcs_median_counts': (CountMetric, {'kwargs': {'always_active': True}, 'prefixes': ['references', 'subsample_types', 'subsample_depths']}),


}

class Reporter:
    def __init__(self, umi_length=None, primers=None,
                 reference_path=None, high_conf_mapq=None, chroms=None,
                 barcode_whitelist=None, barcode_summary=None,
                 gem_groups=None,
                 metrics_dict=METRICS,
                 gene_index=None,
                 umi_min_qual_threshold=0,
                 subsample_types=None,
                 subsample_depths=None,
                 genomes=None):
        self.chroms = chroms
        self.umi_length = umi_length
        self.high_conf_mapq = high_conf_mapq
        self.barcode_whitelist = barcode_whitelist
        self.barcode_summary = barcode_summary
        self.gem_groups = gem_groups
        self.gene_index = gene_index
        self.umi_min_qual_threshold = umi_min_qual_threshold
        random.seed(0)

        self.barcode_types = cr_constants.BARCODE_TYPES
        self.molecule_types = cr_constants.MOLECULE_TYPES
        self.read_types = cr_constants.READ_TYPES
        self.dupe_types = cr_constants.DUPE_TYPES
        self.umi_properties = cr_constants.UMI_PROPERTIES
        self.barcode_properties = cr_constants.BARCODE_PROPERTIES

        self.regions = cr_constants.REGIONS

        assert not (reference_path is not None and genomes is not None)

        if reference_path is not None:
            self.genomes = cr_utils.get_reference_genomes(reference_path)
            self.references = self.genomes + [cr_constants.MULTI_REFS_PREFIX]
            self.has_multiple_genomes = len(self.genomes) > 1
        elif genomes is not None:
            self.genomes = genomes
            self.references = self.genomes + [cr_constants.MULTI_REFS_PREFIX]
            self.has_multiple_genomes = len(self.genomes) > 1
        else:
            self.genomes = []
            self.references = []
            self.has_multiple_genomes = False


        self.subsample_types = subsample_types
        self.subsample_depths = subsample_depths

        if primers is not None:
            self.primers = {}
            self.primer_names = [primer.name for primer in primers]
            for primer_name, primer in zip(self.primer_names, primers):
                seq = primer.seq
                seq_rc = tk_seq.get_rev_comp(primer.seq)
                self.primers[primer_name] = {
                    'seq': seq,
                    'seq_rc': seq_rc,
                }
        else:
            self.primers = {}
            self.primer_names = []

        if umi_length is not None:
            self.umi_primers = {}
            self.poly_t = "T" * min(self.umi_length, cr_constants.UMI_POLYT_SUFFIX_LENGTH)
            for name, seq in self.primers.iteritems():
                self.umi_primers[name] = {
                    'seq': seq['seq'][:self.umi_length],
                    'seq_rc': seq['seq_rc'][:self.umi_length],
                }
        else:
            self.umi_primers = {}
            self.poly_t = None

        self.homopolymers = {}
        self.nucs = tk_seq.NUCS
        for nuc in self.nucs:
            self.homopolymers[nuc] = nuc * cr_constants.HOMOPOLYMER_LENGTH

        self.metrics_dict = metrics_dict
        self.metrics_data = {}
        prefixes = []
        for base, (metric_cls, metric_dict) in self.metrics_dict.iteritems():
            # Track prefixes for passing to the webshim
            prefixes += metric_dict.get('prefixes', [])
            args = metric_dict.get('args', [])
            kwargs = metric_dict.get('kwargs', {})
            for key in self._get_metric_keys(base):
                self.metrics_data[key] = metric_cls(*args, **kwargs)

        self.prefixes = {}
        for prefix in prefixes:
            self.prefixes[prefix] = self.__dict__.get(prefix, [])

        self.reference_path = reference_path
        self.metadata = {}


    def _get_metric_keys(self, name):
        metric_cls, metric_dict = self.metrics_dict[name]
        prefixes = metric_dict.get('prefixes', [])
        kwargs = metric_dict.get('kwargs', {})

        always_active = kwargs.get('always_active', False)

        parts = [[name]]
        for prefix in prefixes:
            prefix = getattr(self, prefix)

            if prefix:
                parts.append(prefix)

        # Check to make sure all specified metrics are present for metrics that are always active
        if always_active and len(parts) != len(prefixes) + 1:
            return []

        # Return the set of keys
        keys = set(itertools.product(*parts))

        # Add bare keys
        keys.add((name,))

        return keys

    def _get_metric_attr(self, *args):
        return self.metrics_data[tuple(args)]

    def get_all_prefixes(self):
        return self.prefixes

    def _set_seq_qual_metrics(self, seq, qual, bases_with_q30_frac, N_bases_frac):
        qvs = tk_fasta.get_qvs(qual)

        num_bases_q30 = np.count_nonzero(qvs >= 30)
        # Don't count no-calls towards Q30 denominator.
        # Assume no-calls get Q <= 2
        num_bases_called = np.count_nonzero(qvs > 2)
        bases_with_q30_frac.add(num_bases_q30, filter=True)
        bases_with_q30_frac.add(num_bases_called-num_bases_q30, filter=False)

        count = seq.count('N')
        N_bases_frac.add(count, filter=True)
        N_bases_frac.add(len(seq)-count, filter=False)

    def raw_fastq_cb(self, read1, read2, bc_read, si_read, umi_read, gem_group):
        get_metric = self._get_metric_attr

        name, seq, qual = read1
        _, seq2, qual2 = read2
        _, bc_seq, bc_qual = bc_read
        _, si_seq, si_qual = si_read
        _, umi_seq, umi_qual = umi_read

        self._set_seq_qual_metrics(seq, qual,
                                   get_metric('read_bases_with_q30_frac'),
                                   get_metric('read_N_bases_frac'))
        get_metric('total_reads').add(1)
        get_metric('total_reads_per_gem_group', gem_group).add(1)

        if seq2:
            self._set_seq_qual_metrics(seq2, qual2,
                                       get_metric('read2_bases_with_q30_frac'),
                                       get_metric('read2_N_bases_frac'))
            get_metric('total_read_pairs').add(1)

        if bc_seq:
            self._set_seq_qual_metrics(bc_seq, bc_qual,
                           get_metric('bc_bases_with_q30_frac'),
                           get_metric('bc_N_bases_frac'))

        if si_seq:
            self._set_seq_qual_metrics(si_seq, si_qual,
                                       get_metric('sample_index_bases_with_q30_frac'),
                                       get_metric('sample_index_N_bases_frac'))
        if umi_seq:
            self._set_seq_qual_metrics(umi_seq, umi_qual,
                           get_metric('umi_bases_with_q30_frac'),
                           get_metric('umi_N_bases_frac'))

        seq_prefix = seq[:cr_constants.READ_PREFIX_LENGTH]
        get_metric('top_read_prefixes').add(seq_prefix)
        if seq2:
            seq2_prefix = seq2[:cr_constants.READ_PREFIX_LENGTH]
            get_metric('top_read2_prefixes').add(seq2_prefix)

        read_with_perfect_primer_or_homopolymer = False
        for nuc, homopolymer_seq in self.homopolymers.iteritems():
            perfect_homopolymer = homopolymer_seq in seq
            read_with_perfect_primer_or_homopolymer = (
                read_with_perfect_primer_or_homopolymer or perfect_homopolymer)

            perfect_homopolymer_frac = get_metric('perfect_homopolymer_frac', nuc)
            perfect_homopolymer_frac.add(1, filter=perfect_homopolymer)

        for primer_name, primer in self.primers.iteritems():
            primer_seq = primer['seq']
            primer_seq_rc = primer['seq_rc']
            present = primer_seq in seq or primer_seq_rc in seq

            perfect_primer = present > 0
            read_with_perfect_primer_or_homopolymer = (
                read_with_perfect_primer_or_homopolymer or perfect_primer)

            perfect_primer_frac = get_metric('perfect_primer_frac', primer_name)
            perfect_primer_frac.add(1, filter=present)

        get_metric('read_with_perfect_primer_or_homopolymer_frac').add(
            1, filter=read_with_perfect_primer_or_homopolymer)

    def _set_mapping_metrics(self, read_type, read_type_set):
        for genome, region in itertools.product(self.genomes, self.regions):
            is_read_type = (genome, region) in read_type_set
            reads_frac = self._get_metric_attr('reads_frac', genome, region, read_type)
            reads_frac.add(1, filter=is_read_type)

        for region in self.regions:
            is_read_type = any([(genome, region) in read_type_set for genome in self.genomes])
            multi_reads_frac = self._get_metric_attr('reads_frac', cr_constants.MULTI_REFS_PREFIX, region,
                                                     read_type)
            multi_reads_frac.add(1, filter=is_read_type)

    def _set_antisense_metrics(self, genome_set):
        for genome in self.genomes:
            antisense_reads_frac = self._get_metric_attr('antisense_reads_frac', genome)
            antisense_reads_frac.add(1, filter=genome in genome_set)
        multi_antisense_reads_frac = self._get_metric_attr('antisense_reads_frac', cr_constants.MULTI_REFS_PREFIX)
        multi_antisense_reads_frac.add(1, filter=len(genome_set) > 0)

    def _set_barcode_reads_metrics(self, read_type, read_type_set, bc):
        for genome in self.genomes:
            is_read_type = (genome, cr_constants.TRANSCRIPTOME_REGION) in read_type_set
            if is_read_type:
                barcode_reads = self._get_metric_attr(
                    'barcode_reads', genome, cr_constants.TRANSCRIPTOME_REGION, read_type)
                barcode_reads.add(bc)

        # Don't always-report the multi prefix for the barcode_reads metrics
        if self.has_multiple_genomes:
            is_read_type = any([(genome, cr_constants.TRANSCRIPTOME_REGION) in read_type_set for genome in self.genomes])
            if is_read_type:
                multi_barcode_reads = self._get_metric_attr(
                    'barcode_reads', cr_constants.MULTI_REFS_PREFIX,
                    cr_constants.TRANSCRIPTOME_REGION, read_type)
                multi_barcode_reads.add(bc)

    def raw_umi_cb(self, seq, qual):
        """ Returns the processed sequence """
        get_metric = self._get_metric_attr
        property_metric = 'umi_property_frac'
        filter_metric = 'umi_filter_frac'

        low_min_qual = cr_utils.min_qual_below(qual, self.umi_min_qual_threshold)
        has_n = 'N' in seq
        is_homopolymer = cr_utils.is_homopolymer_seq(seq)
        has_primer = cr_utils.find_any_primers(seq, self.umi_primers)
        has_polyt = seq[:-cr_constants.UMI_POLYT_SUFFIX_LENGTH] == self.poly_t

        get_metric(property_metric, cr_constants.LOW_MIN_QUAL_UMI_FILTER).add(1, low_min_qual)
        get_metric(property_metric, cr_constants.HAS_N_UMI_FILTER).add(1, has_n)
        get_metric(property_metric, cr_constants.HOMOPOLYMER_UMI_FILTER).add(1, is_homopolymer)
        get_metric(property_metric, cr_constants.PRIMER_UMI_FILTER).add(1, has_primer)
        get_metric(property_metric, cr_constants.POLYT_UMI_FILTER).add(1, has_polyt)

        valid = True

        get_metric(filter_metric, cr_constants.LOW_MIN_QUAL_UMI_FILTER).add(1, low_min_qual)
        valid = valid and not low_min_qual

        get_metric(filter_metric, cr_constants.HAS_N_UMI_FILTER).add(1, valid and has_n)
        valid = valid and not has_n

        get_metric(filter_metric, cr_constants.HOMOPOLYMER_UMI_FILTER).add(1, valid and is_homopolymer)
        valid = valid and not is_homopolymer

        return seq if valid else None

    def raw_barcode_cb(self, seq, qual, whitelist):
        """ Returns the processed sequence """
        get_metric = self._get_metric_attr
        property_metric = 'barcode_property_frac'
        filter_metric = 'barcode_filter_frac'

        missed_whitelist = not cr_utils.is_barcode_on_whitelist(seq, whitelist)
        has_n = 'N' in seq
        is_homopolymer = cr_utils.is_homopolymer_seq(seq)
        low_min_qual = cr_utils.min_qual_below(qual, cr_constants.BARCODE_MIN_QUAL_THRESHOLD)

        get_metric(property_metric, cr_constants.MISS_WHITELIST_BARCODE_FILTER).add(1, missed_whitelist)
        get_metric(property_metric, cr_constants.HAS_N_BARCODE_FILTER).add(1, has_n)
        get_metric(property_metric, cr_constants.HOMOPOLYMER_BARCODE_FILTER).add(1, is_homopolymer)
        get_metric(property_metric, cr_constants.LOW_MIN_QUAL_BARCODE_FILTER).add(1, low_min_qual)

        valid = True

        get_metric(filter_metric, cr_constants.MISS_WHITELIST_BARCODE_FILTER).add(1, missed_whitelist)
        valid = valid and not missed_whitelist

        return seq if valid else None

    def aligned_bam_cb(self, genome_alignments, transcriptome_alignments,
                       barcode_info=None, umi_info=None, any_antisense=False):
        assert self.high_conf_mapq
        get_metric = self._get_metric_attr

        read_mapped, read_conf_mapped, read_antisense = set(), set(), set()
        for alignment in genome_alignments:
            genome = cr_utils.get_genome_from_read(alignment, self.chroms, self.genomes)

            for region in self.regions:
                is_mapped = cr_utils.is_read_mapped(alignment, region, gene_index=self.gene_index,
                                                    chroms=self.chroms)
                if is_mapped:
                    read_mapped.add((genome, region))

                is_conf_mapped = cr_utils.is_read_conf_mapped(alignment, region,
                                                              self.high_conf_mapq,
                                                              gene_index=self.gene_index,
                                                              chroms=self.chroms)
                if is_conf_mapped:
                    read_conf_mapped.add((genome, region))
                    if region == cr_constants.TRANSCRIPTOME_REGION:
                        if any_antisense:
                            read_antisense.add(genome)

        unmapped_metric = self._get_metric_attr('unmapped_reads_frac')
        unmapped_metric.add(1, filter=len(read_mapped) == 0)

        self._set_mapping_metrics(cr_constants.MAPPED_READ_TYPE, read_mapped)
        self._set_mapping_metrics(cr_constants.CONF_MAPPED_READ_TYPE, read_conf_mapped)
        self._set_antisense_metrics(read_antisense)

        # Compute insert size
        if cr_utils.is_read_conf_mapped(genome_alignments[0], cr_constants.TRANSCRIPTOME_REGION,
                                        self.high_conf_mapq,
                                        gene_index=self.gene_index, chroms=self.chroms):
            if any(alignment.is_read1 or alignment.is_read2 for alignment in genome_alignments):
                # Compute paired end insert size
                insert_size = cr_stats.compute_pe_insert_size(genome_alignments[0],
                                                              transcriptome_alignments)
            else:
                # Infer 3' single end insert size
                insert_size = cr_stats.compute_se_3p_insert_size(genome_alignments[0],
                                                                 self.gene_index)
            self._set_insert_size_metrics(genome_alignments[0], insert_size)

        # Barcode and umi metrics
        if umi_info:
            get_metric('good_umi_frac').add(1, filter=umi_info.processed_seq)
            get_metric('top_raw_umis').add(umi_info.raw_seq)

            if umi_info.processed_seq:
                get_metric('effective_umi_diversity').add(umi_info.processed_seq)
                get_metric('top_processed_umis').add(umi_info.processed_seq)

            read_conf_mapped_umi = read_conf_mapped and umi_info.processed_seq
            if read_conf_mapped_umi:
                genome, _ = next(iter(read_conf_mapped))
                for reference in [genome, cr_constants.MULTI_REFS_PREFIX]:
                    conf_mapped_effective_umi_diversity = get_metric(
                        'effective_umi_diversity', reference, cr_constants.CONF_MAPPED_READ_TYPE)
                    conf_mapped_effective_umi_diversity.add(umi_info.processed_seq)

                    conf_mapped_top_processed_umis = get_metric(
                        'top_processed_umis', reference, cr_constants.CONF_MAPPED_READ_TYPE)
                    conf_mapped_top_processed_umis.add(umi_info.processed_seq)

        if barcode_info:
            get_metric('good_bc_frac').add(1, filter=barcode_info.processed_seq)
            get_metric('corrected_bc_frac').add(1, filter=cr_utils.is_barcode_corrected(barcode_info.raw_seq,
                                                                                        barcode_info.processed_seq))
            get_metric('top_raw_barcodes').add(barcode_info.raw_seq)
            if barcode_info.processed_seq:
                get_metric('barcodes_detected').add(barcode_info.processed_seq)

            read_conf_mapped_bc = read_conf_mapped if read_conf_mapped and barcode_info.processed_seq else set()
            self._set_mapping_metrics(cr_constants.CONF_MAPPED_BC_READ_TYPE, read_conf_mapped_bc)

            if barcode_info.processed_seq:
                get_metric('effective_barcode_diversity').add(barcode_info.processed_seq)
                get_metric('barcode_reads').add(barcode_info.processed_seq)

                self._set_barcode_reads_metrics(cr_constants.MAPPED_READ_TYPE, read_mapped, barcode_info.processed_seq)
                self._set_barcode_reads_metrics(cr_constants.CONF_MAPPED_READ_TYPE, read_conf_mapped, barcode_info.processed_seq)
                self._set_barcode_reads_metrics(cr_constants.CONF_MAPPED_BC_READ_TYPE, read_conf_mapped_bc, barcode_info.processed_seq)

            if read_conf_mapped_bc:
                genome, _ = next(iter(read_conf_mapped))
                for reference in [genome, cr_constants.MULTI_REFS_PREFIX]:
                    conf_mapped_effective_barcode_diversity = get_metric(
                        'effective_barcode_diversity', reference, cr_constants.CONF_MAPPED_READ_TYPE)
                    conf_mapped_effective_barcode_diversity.add(barcode_info.processed_seq)

    def transcript_primer_metrics_cb(self, genome_reads, transcriptome_reads, transcripts, gene_index):
        """
        Calculates primer position metrics within transcriptome aligned reads.
        # TODO: Move to PD
        # TODO eventually would like to support primer matches and paired-end calculation

        Args:
            transcriptome_reads (list of pysam.AlignedRead): a list of transcriptome aligned reads
            genome_reads (list of pysam.AlignedRead): a list of genome aligned reads
            transcripts (dict): Map of BAM transcript IDs to trancript IDs compatible with gene_index
            gene_index (GeneIndex): gene index object built from GTF file

        """

        distances_from_three_prime = []
        distances_from_five_prime = []

        # Check if the genome read is confidently aligned
        confidently_mapped = cr_utils.is_read_conf_mapped(genome_reads[0], cr_constants.TRANSCRIPTOME_REGION, self.high_conf_mapq)

        # Track longest transcript(s) and only calculate final metric using this transcript
        longest_transcript_length = 0

        for transcriptome_read in transcriptome_reads:

            # Only calculate transcript position metrics on confidently mapped reads
            if confidently_mapped:
                # Note: for paired-end reads, only consider read1
                assert len(genome_reads) == 1 \
                    or len(genome_reads) == 2 and genome_reads[0].is_read1 \
                    or genome_reads[0].has_tag(cr_constants.MULTIMAPPER_TAG)
                transcript = gene_index.transcripts[transcripts[transcriptome_read.tid]]
                transcript_length = transcript.length

                # Reset tracking if find a longer compatible transcript
                if transcript_length > longest_transcript_length:
                    distances_from_three_prime = []
                    distances_from_five_prime = []

                # Add distances to tracking for longest transcript
                if transcript_length >= longest_transcript_length:
                    longest_transcript_length = transcript_length

                    distance_to_start = transcriptome_read.pos  # position relative to transcript start
                    distance_to_end = transcript_length - (transcriptome_read.pos + transcriptome_read.qlen)

                    distances_from_three_prime.append(max(distance_to_end, 0))
                    distances_from_five_prime.append(max(distance_to_start, 0))

        # Take median of all positions derived from the longest transcript(s)
        if len(distances_from_three_prime) > 0:
            three_prime_distance = self._get_metric_attr('read_distance_to_transcript_three_prime')
            three_prime_distance.add(np.median(np.array(distances_from_three_prime)))

            five_prime_distance = self._get_metric_attr('read_distance_to_transcript_five_prime')
            five_prime_distance.add(np.median(np.array(distances_from_five_prime)))


    def _set_dupe_metrics(self, read, dupe_type, genome):
        dupe_reads_frac = self._get_metric_attr('dupe_reads_frac', genome, dupe_type)
        dupe_reads_rate = self._get_metric_attr('dupe_reads_rate', genome, dupe_type)
        dupe_reads_frac.add(1, filter=read.is_duplicate)
        dupe_reads_rate.add(1, numerator=read.is_duplicate)

    def mark_dupes_bam_cb(self, read, dupe_type):
        assert self.high_conf_mapq
        if not cr_utils.is_read_dupe_candidate(read, self.high_conf_mapq):
            return

        genome = cr_utils.get_genome_from_read(read, self.chroms, self.genomes)
        self._set_dupe_metrics(read, dupe_type, genome)
        self._set_dupe_metrics(read, dupe_type, cr_constants.MULTI_REFS_PREFIX)

    def mark_dupes_group_cb(self, gene_id, umis, dupe_type):
        total_counts = sum(umis.values())
        total_umis = len(umis)
        if any([count > 1 for count in umis.itervalues()]):
            umi_hamming_distance = 0
        else:
            umi_hamming_distance = cr_utils.get_kmers_hamming_distance(umis.keys())

        for reference in [cr_utils.get_genome_from_str(gene_id, self.genomes), cr_constants.MULTI_REFS_PREFIX]:
            if total_counts > 0:
                reads_per_dupe_group_histogram = self._get_metric_attr(
                    'reads_per_dupe_group_histogram', reference, dupe_type)
                reads_per_dupe_group_histogram.add(total_counts)

            if total_umis > 0:
                umis_per_dupe_group_histogram = self._get_metric_attr(
                    'umis_per_dupe_group_histogram', reference, dupe_type)
                umis_per_dupe_group_histogram.add(total_umis)

            reads_per_molecule_histogram = self._get_metric_attr(
                'reads_per_molecule_histogram', reference, dupe_type)
            for count in umis.itervalues():
                reads_per_molecule_histogram.add(count)

            if umi_hamming_distance is not None:
                umi_hamming_distance_per_dupe_group_histogram = self._get_metric_attr(
                    'umi_hamming_distance_per_dupe_group_histogram', reference, dupe_type)
                umi_hamming_distance_per_dupe_group_histogram.add(umi_hamming_distance)

    def mark_dupes_corrected_cb(self, read):
        if read.is_secondary:
            return

        raw_umi_seq = cr_utils.get_read_raw_umi(read)
        processed_umi_seq = cr_utils.get_read_umi(read)
        self._get_metric_attr('corrected_umi_frac').add(1, filter=cr_utils.is_umi_corrected(raw_umi_seq, processed_umi_seq))

    def count_genes_bam_cb(self, read, use_umis=True):
        assert self.high_conf_mapq

        if read.is_secondary:
            return False, None, None, None

        if use_umis:
            conf_mapped_deduped = cr_utils.is_read_conf_mapped_to_transcriptome_deduped(read, self.high_conf_mapq)
        else:
            conf_mapped_deduped = cr_utils.is_read_conf_mapped_to_transcriptome_barcoded(read, self.high_conf_mapq)

        genome = cr_utils.get_genome_from_read(read, self.chroms, self.genomes)
        for reference in self.references:
            conf_mapped_deduped_frac = self._get_metric_attr(
                'reads_frac', reference,
                cr_constants.TRANSCRIPTOME_REGION,
                cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE)
            conf_mapped_deduped_frac.add(1, filter=conf_mapped_deduped and reference in [genome, cr_constants.MULTI_REFS_PREFIX])

        if conf_mapped_deduped:
            bc = cr_utils.get_read_barcode(read)
            umi = cr_utils.get_read_umi(read)
            for reference in [genome, cr_constants.MULTI_REFS_PREFIX]:
                if bc:
                    conf_mapped_deduped_effective_barcode_diversity = self._get_metric_attr(
                        'effective_barcode_diversity', reference,
                        cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE)
                    conf_mapped_deduped_effective_barcode_diversity.add(bc)

                    # Only report barcode_reads for multi_* if there are multiple genomes
                    if reference != cr_constants.MULTI_REFS_PREFIX or self.has_multiple_genomes:
                        conf_mapped_deduped_barcode_reads = self._get_metric_attr(
                            'barcode_reads', reference,
                            cr_constants.TRANSCRIPTOME_REGION,
                            cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE)
                        conf_mapped_deduped_barcode_reads.add(bc)

                if umi:
                    conf_mapped_deduped_effective_umi_diversity = self._get_metric_attr(
                        'effective_umi_diversity', reference,
                        cr_constants.CONF_MAPPED_DEDUPED_READ_TYPE)
                    conf_mapped_deduped_effective_umi_diversity.add(umi)

            gene_ids = cr_utils.get_read_gene_ids(read)
            assert len(gene_ids) == 1
            gene_id = gene_ids[0]

            return True, genome, gene_id, bc

        return False, None, None, None

    def _set_insert_size_metrics(self, read, insert_size):
        genome = cr_utils.get_genome_from_read(read, self.chroms, self.genomes)

        if insert_size is not None:
            for reference in [genome, cr_constants.MULTI_REFS_PREFIX]:
                median_insert_size = self._get_metric_attr('median_insert_size', reference)
                insert_size_histogram = self._get_metric_attr('insert_size_histogram', reference)
                median_insert_size.add(insert_size)
                insert_size_histogram.add(insert_size)

                iqr_insert_size = self._get_metric_attr('iqr_insert_size', reference)
                iqr_insert_size.add(insert_size)

    def subsampled_duplication_frac_cb(self, subsampled_raw_matrices, mol_counter,
                                       subsampled_rate, subsampled_type, subsampled_depth,
                                       mapped_reads):
        """
        Calculates the subsampled duplication rate for a subsampled raw matrix.
        NOTE: This assumes the sampled depth is identical to the targeted depth.
              It will be inaccurate at low molecule counts.

        Args:
            subsampled_raw_matrices (GeneBCMatrices): subsampled GeneBCMatrices object (must be raw matrix)
            mol_counter (MoleculeCounter): molecule counter object for the subsampled matrix
            subsampled_rate (float): rate of subsampling used
            subsampled_type (str): type of subsampling performed
            subsampled_depth (int): target reads per cell for subsampled matrix
            mapped_reads (int): total mapped reads after subsampling

        """

        for reference in self.references:
            subsampled_duplication_frac = self._get_metric_attr('subsampled_duplication_frac', reference, subsampled_type, subsampled_depth)

            if reference == cr_constants.MULTI_REFS_PREFIX:
                total_molecules = subsampled_raw_matrices.get_reads_per_bc().sum()
            else:
                total_molecules = subsampled_raw_matrices.matrices[reference].m.sum()

            subsampled_duplication_frac.set_value(numerator=mapped_reads - total_molecules,
                                                  denominator=mapped_reads)

    def summarize_subsampled_matrices_cb(self, filtered_mats, subsample_type, subsample_depth):
        """
        Computes simple summary metrics such as median genes detected and UMI counts on subsampled filtered matrices

        Args:
            filtered_mats (GeneBCMatrices): subsampled and filtered GeneBCMatrices
            subsample_type (string): subsampling type
            subsample_depth (int): target depth per cell for subsampling

        """
        for genome in self.genomes:
            if filtered_mats is not None:
                matrix = filtered_mats.matrices[genome]
                genes_detected = np.median(matrix._sum(matrix.m >= cr_constants.MIN_READS_PER_GENE, axis=0))
                median_counts = np.median(matrix._sum(matrix.m, axis=0))

            subsampled_filtered_bc_median_unique_genes_detected = self._get_metric_attr('subsampled_filtered_bcs_median_unique_genes_detected', genome, subsample_type, subsample_depth)
            subsampled_filtered_bc_median_unique_genes_detected.set_value(genes_detected)

            subsampled_filtered_bcs_median_counts = self._get_metric_attr('subsampled_filtered_bcs_median_counts', genome, subsample_type, subsample_depth)
            subsampled_filtered_bcs_median_counts.set_value(median_counts)



    def store_pipeline_metadata(self):
        if self.metadata is None:
            self.metadata = {}
        self.metadata[cr_constants.CELLRANGER_VERSION_KEY] = cr_utils.get_version()

    def store_reference_metadata(self):
        if self.metadata is None:
            self.metadata = {}
        if not self.reference_path:
            return

        ref_metadata = cr_utils._load_reference_metadata_file(self.reference_path)
        for key in cr_constants.REFERENCE_METADATA_KEYS:
            value = ref_metadata.get(key, '')
            if value is None:
                value = ''
            if np.isscalar(value):
                self.metadata['reference_%s' % key] = value
            else:
                self.metadata['reference_%s' % key] = ', '.join(str(x) for x in value)

    def store_chemistry_metadata(self, chemistry_def):
        """ Store the chemistry definition as metrics in the summary json """
        if self.metadata is None:
            self.metadata = {}
        for key, value in chemistry_def.iteritems():
            self.metadata['chemistry_%s' % key] = value

    def merge(self, reporter):
        for key, metric2 in reporter.metrics_data.iteritems():
            metric1 = self.metrics_data.get(key)
            if isinstance(metric2, Metric):
                if metric1 and metric1.active and metric2.active:
                    metric1.merge(metric2)
                elif metric2.active:
                    self.metrics_data[key] = metric2

    def report(self, report_type):
        summary = {}
        for key, metric in self.metrics_data.iteritems():
            if metric.active and metric.report_type == report_type:

                name = '_'.join([str(x) for x in key[1:] + (key[0],)]) # ensure keys are strings before output
                summary[name] = metric.report()

        if self.metadata is not None and report_type == cr_constants.DEFAULT_REPORT_TYPE:
            summary.update(self.metadata)

        return summary

    def report_summary_json(self, filename):
        data = self.report(cr_constants.DEFAULT_REPORT_TYPE)
        with open(filename, 'w') as f:
            tk_safe_json.dump_numpy(tk_safe_json.json_sanitize(data), f, pretty=True)

    def report_barcodes_h5(self, filename):
        data = self.report('barcodes')

        if self.barcode_whitelist:
            bc_sequences = cr_utils.format_barcode_seqs(self.barcode_whitelist, self.gem_groups)
        elif self.barcode_summary is not None:
            bc_sequences = self.barcode_summary
        else:
            # Get all observed bc sequences
            bc_sequences = sorted(list(reduce(lambda x, y: x.union(y.keys()), data.values(), set())))

        # Build the columns for the table
        bc_table_cols = {cr_constants.H5_BC_SEQUENCE_COL: bc_sequences}

        for metric_name, metric_data in data.iteritems():
            counts = np.array([metric_data.get(bc, 0) for bc in bc_sequences], dtype=np.uint32)
            bc_table_cols[metric_name] = counts

        cr_utils.write_h5(filename, bc_table_cols)


    def save(self, filename):
        with open(filename, 'wb') as f:
            cPickle.dump(self, f, protocol=cPickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(filename):
        with open(filename, 'rb') as f:
            return cPickle.load(f)

def merge_reporters(filenames):
    reporter = None
    for filename in filenames:
        tmp_reporter = Reporter.load(filename)
        if reporter is None:
            reporter = tmp_reporter
        else:
            reporter.merge(tmp_reporter)

    return reporter

def merge_jsons(in_filenames, out_filename, dicts=[]):
    """ Merge a list of json files and optional dicts """
    d = cr_utils.merge_jsons_as_dict(in_filenames)
    for data in dicts:
        cr_utils.update_require_unique_key(d, data)

    with open(out_filename, 'w') as f:
        json.dump(tk_safe_json.json_sanitize(d), f, indent=4, sort_keys=True)

def merge_h5(in_filenames, out_filename):
    """ Merge a list of h5 files """
    out_h5 = h5.File(out_filename, 'a')
    for filename in in_filenames:
        if filename is None:
            continue
        in_h5 = h5.File(filename, 'r')
        for name in in_h5.keys():
            # If the dataset already exists,
            # They must be equal or one must be all-zero.
            if name in out_h5.keys():
                src_data, dest_data = in_h5[name][()], out_h5[name][()]
                if src_data.dtype.kind != 'S' and dest_data.dtype.kind != 'S':
                    # Both numeric
                    if not np.any(src_data):
                        # Source is all zero. Do nothing.
                        continue
                    elif not np.any(dest_data):
                        # Dest is all zero. Overwrite.
                        del out_h5[name]
                        h5.h5o.copy(in_h5.id, name, out_h5.id, name)
                    else:
                        # Both non-zero. Assert equality and do nothing.
                        assert np.array_equal(src_data, dest_data)
                else:
                    # Either are non-numeric. Assert equality and do nothing.
                    assert np.array_equal(src_data, dest_data)
            else:
                # Only exists in src. Copy to dest.
                h5.h5o.copy(in_h5.id, name, out_h5.id, name)

    out_h5.flush()
    out_h5.close()
