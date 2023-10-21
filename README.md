# EBFilter (Empirical Bayesian Mutation Filtering)

## Introduction

EBFilter is a software for filtering false positives somatic mutations in cancer genome sequencing data analysis.
EBFilter accepts a list of candidates somatic mutation already narrowed down to some extents, and performs the main step of [EBCall](https://github.com/friend1ws/EBCall), 

1. estimate the parameters of the beta-binomial sequencing error model using multiple non-matched control sequencing data at the position of interest
2. get the predictive mismatch ratio derived from the above estimation and compare it with the observed mismatch ratio of tumor samples
3. if the mismatch ratio of the tumor sample is significantly deviated from the predicted mismatch ratio, then we identify it as highly-likely somatic mutation

Therefore, you can use EBFilter along with your own mutation calling program  
or some popular mutation callers (e.g., MuTect, VarScan 2),
which we believe will reduce large parts of false positives with a slight computational cost.

## Paper

We would like you to kindly cite the following paper when you use this software; 

"[An empirical Bayesian framework for mutation detection from cancer genome sequencing data](http://nar.oxfordjournals.org/content/41/7/e89.long)", Shiraishi et al.,  
Nucleic Acids Research, 2013.


# Motivation

One major source of false positive somatic mutation in cancer genome sequencing data analysis is,
"broken symmetry of mismatches between tumor and normal sequencing data that randomly occurs at the error-prone sites"
(see e.g., Figure 1, Shiraishi et al., NAR, 2013).
Therefore, one of the most important techniques for reducing false positive mutations is to check many control sequencing data (including non-matched normal samples), to see whether the candidates are artifacts produced at the error-prone sites, and many somatic mutations calling pipeline adopt this strategy.

EBCall uses 

1. estimate the parameters of the sequencing error model using multiple non-matched contorol sequencing data at the position of interest
2. get the predictive mismatch ratio obtained from the above estimation and compare it with the observed mismatch ratio of tumor samples
3. if the mismatch ratio of the tumor sample is significantly deviated from the predicted mismatch ratio, then we idenitfy it as highly-likely somatic mutation


EBCall is actually integrating many other steps (Fisher's exact test), and could not perform purely the above beta-binomial filteringt step.
Therefore, we decided to implement a software which just perform beta-binomial step for already collected candidates of somatic mutations.

Therefore, you can use this software after performing popular mutation callers (e.g., mutects, VarScan2,and so on),
or your own inhouse mutation calling program, which we believe will reduce large parts of false positives.

# Features of this version.
This version is modified for use in ogawa-lab and includes the following changes.
- Stable results are produced even with high depth samples.
- Optimised handling of overlapping reads.
- Reduced file access load and increased speed.


## Dependency

### Software
[samtools](http://www.htslib.org/)

### Python
Python (>= 2.7), `pysam`, `scipy`, `numpy` packages


## Preparation
- add path to samtools.
- **target somatic mutation candidats**: the somatic mutation candidates (should be .vcf format).
- **target tumor sample**: the indexed bam file of the target tumor sample.
- **list of normal reference samples**: the list of paths to the indexed bam files for non-paired normal reference samples. Please name the text file as you like (e.g., myNormalRef.txt), and list the paths of .bam files as follows:  

		/home/yshira/ngs/data/sequence/normalreference1.bam
		/home/yshira/ngs/data/sequence/normalreference2.bam
		...
		/home/yshira/ngs/data/sequence/normalreference10.bam

## Commands
    EBFilter_0.2.5 [-h] [--version] [-t thread_num]
                [-q mapping_qual_thres] [-Q base_qual_thres] [-ff filter_flags]
                [--loption] [--region REGION] [--epsilon EPSILON] [--method METHOD] [--debug]
                target.tsv target.bam controlBam_list.txt output.tsv





