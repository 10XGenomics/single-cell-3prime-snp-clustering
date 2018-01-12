# Single Cell 3' RNA SNP Clustering Pipeline (`snpclust`)

Using single cell RNA-sequencing data produced by the 10x Genomics Chromium Single Cell 3' Solution, cluster cells into two groups by their single nucleotide polymorphisms (SNPs). The pipeline takes as input a BAM file and barcode list produced by Cell Ranger and produces cluster calls for each cell.

**While this software was used to analyze the donor-host single cell RNA-seq data in AML patients, this software has not been extensively tested on other sample types or applications. This software is not supported by 10x Genomics.**

**This software is for research use only. The code in this repository is licensed under [AGPLv2](http://www.affero.org/agpl2.html).**

## System Requirements
Same as [Cell Ranger system requirements](https://support.10xgenomics.com/single-cell/software/overview/system-requirements)

## Dependencies
- [Ranger (v1.0.1 ONLY - newer versions may not work)](https://support.10xgenomics.com/developers/software)
- [Cell Ranger reference data package](https://support.10xgenomics.com/single-cell/software/downloads/latest)

## Install
1. Install Ranger
  * The example assumes you have untar'ed Ranger to /opt/10x/ranger-1.0.1

2. Install the cellranger reference data
  - Install the appropriate Cell Ranger reference data.
    * The example assumes you have installed the hg19-1.2.0 refdata to /opt/10x/refdata-cellranger-hg19-1.2.0

3. Clone this repository to. e.g., /opt/10x/cellranger-snpclust-src-1.0.0


# Run the pipeline

## Download example data
```
cd /opt/10x
# BAM file (WARNING: 9GB file)
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/frozen_pbmc_b_c_90_10/frozen_pbmc_b_c_90_10_possorted_genome_bam.bam -O example.bam

# Get list of cell barcodes
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/frozen_pbmc_b_c_90_10/frozen_pbmc_b_c_90_10_filtered_gene_bc_matrices.tar.gz -O matrices.tar.gz
tar -zxf matrices.tar.gz
cp filtered_matrices_mex/hg19/barcodes.tsv example_barcodes.tsv
```

## Construct an mro file
This software uses the 10x Genomics Martian platform. The pipeline is invoked by supplying a Martian Runtime Object (.mro) file containing input parameters to the Martian Pipeline Runner (`mrp`).

The following arguments are paths to your inputs:
- `reference_path` (path to 10x refdata)
- `possorted_genome_bam` (the BAM file produced by Cell Ranger)
- `cell_barcodes` (a text file containing cell-barcodes that were called as cells by Cell Ranger; one barcode per line)

/opt/10x/example.mro:
```
@include "snp_clusterer_cs.mro"

call SNP_CLUSTERER_CS(
    reference_path       = "/opt/10x/refdata-cellranger-hg19-1.2.0",
    possorted_genome_bam = "/opt/10x/example.bam",
    cell_barcodes        = "/opt/10x/example_barcodes.tsv",
    align                = {"high_conf_mapq":255},
    min_snp_call_qual    = 0,
    min_bcs_per_snp      = 2,
    min_snp_obs          = 1,
    min_snp_base_qual    = 1,
    base_error_rate      = null,
    seed                 = 0,
    iterations           = 10000,
    burn                 = 1000,
    skip                 = 5,
    cluster_ground_truth = null,
    min_depth            = 0,
    min_posterior        = 0.8,
)
```

## Enter the 10x environment and run the pipeline
```
# Assuming ranger was extracted to /opt/10x/ranger-1.0.1
# and this repo was cloned to /opt/10x/cellranger-snpclust-src-1.0.0
export PATH=/opt/10x/ranger-1.0.1:/opt/10x/ranger-1.0.1/freebayes/v1.0.2:$PATH
export MROPATH=/opt/10x/cellranger-snpclust-src-1.0.0/mro:$MROPATH
export PYTHONPATH=/opt/10x/cellranger-snpclust-src-1.0.0/lib/python/:/opt/10x/ranger-1.0.1/ranger-cs/1.0.1/tenkit/lib/python:$PYTHONPATH
ranger mrp example.mro example
```

# Understand the output

The pipeline produces 5 main output files/directories in the outs folder:

- likelihood_allele_bc_matrices_mex: a directory that contains the likelihood of a cell having one of the 3 genotypes (reference/reference  (1/1), alt/alt (0/0), reference/alt(1/0)) across all SNVs. There are 3 files for each genotype: barcodes.tsv (cell barcode Ids in tsv format), genes.tsv (chromosome location of SNVs in tsv format), matrix.mtx (a sparse matrix containing the likelihood of the genotype for each cell at each SNV)

- likelihood_allele_bc_matrices_h5.h5: the same information as in likelihood_allele_bc_matrices_mex, organized in hdf5 format

- raw_allele_bc_matrices_mex:a directory that contains the raw allele counts of each SNV across all cells (consisting of 2 directories, "alt" and "ref" indicating counts for alternative and reference alleles respectively). There are 3 files for each directory: barcodes.tsv (cell barcode Ids in tsv format), genes.tsv (chromosome location of SNVs in tsv format), matrix.mtx (a sparse matrix containing the likelihood of the genotype for each cell at each SNV)

- raw_allele_bc_matrices_h5.h5: the same information as in raw_allele_bc_matrices_mex, organized in hdf5 format

- summary.json: summary of genotype calls over the K=1 and K=2 models. The most important entries include: "model1_call" (the assignment of each cell under K=1 model with 1 indicating that cell=genotype1 is True, and 0 indicating that cell=genotype1 is False), "model1_probability" (the probability that K=1 model is true for a cell), "model2_call" (the assignment of each cell under K=2 model with 1 indicating that cell=genotype1 is True, and 0 indicating that cell=genotype1 is False, i.e. cell=genotype2 True), and "model2_probability" (the probability that K=2 model is true for a cell). We applied a requirement that 90% of the cells have a posterior probability >75% to select the K=2 model over the K=1 model. Otherwise, K=1 model is selected.

To ascertain the donor and host identity of each cell in a sample, summary.json is first used to decide if K=1 (there is 1 genotype group among cells) or K=2 (there are 2 genotype groups among cells) model is used. This can be done by evaluating whether there is >90% of cells with model2_probability>75%. If yes, K=2 model is selected, and model2_call can be used to ascertain the genotype group of each cell. Otherwise, K=1 model is selected, and model1_call is used to ascertain the genotype group of each cell.

To get more information on the genotype of each SNV, files in likelihood_allele_bc_matrices_mex and raw_allele_bc_matrices_mex directories can be used. For example, the list of all the SNVs identified can be found in genes.tsv in one of the directories of raw_allele_bc_matrices_mex. The alternative allele frequency of each SNV can be calculated by using the allele count matrix in each of the directories of raw_allele_bc_matrices_mex. The likelihood of genotype assignment for each SNV in a cell can be obtained from the matrix.mtx in each of the directories of likelihood_allele_bc_matrices_mex. The genotype that received the highest likehood score is inferred to be the genotype of the SNV in the cell. See the following R scripts for examples.

# Example R scripts
```
  # read in summary json file, and decide which model to use
  json_results<-fromJSON(json_str="outs/summary.json")
  sum(json$example_model2_probability>0.75)
  table(json$example_model2_call)

  # read in reference and alt allele counts
  snps<-read.delim("outs/raw_allele_bc_matrices_mex/example_alt/genes.tsv",as.is=T,header=F)[,1]
  ref<-as.matrix(readMM("outs/raw_allele_bc_matrices_mex/example_ref/matrix.mtx"))
  alt<-as.matrix(readMM(outs/raw_allele_bc_matrices_mex/example_alt/matrix.mtx"))
  ref_rs<-rowSums(ref)
  alt_rs<-rowSums(alt)
  snps_count<-rowSums(ref+alt)
  snps_freq<-alt/(alt+ref+1)

  # read in genotype likelihood at each snv in each cell
  rr<-readMM("outs/likelihood_allele_bc_matrices_mex/example_1_0|0/matrix.mtx")
  ra<-readMM(outs/likelihood_allele_bc_matrices_mex/example_1_1|0/matrix.mtx")
  aa<-readMM(outs/likelihood_allele_bc_matrices_mex/example_1_1|1/matrix.mtx"))

  # calculate the inferred genotype at each snv in each cell
  rr_rs<-rowSums(rr)
  ra_rs<-rowSums(ra)
  aa_rs<-rowSums(aa)
  gt<-apply(cbind(rr_rs,ra_rs,aa_rs),1,which.max)

  ```
