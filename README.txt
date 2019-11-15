Dependencies
============

Development and testing was done on an up-to-date Ubuntu Linux 19.10 instance.

Requirements:
- R 3.6.1
- bedtools 2.27.1
- xargs (e.g. from findutils)
- gzip
- grep, sed
- awk (tested with mawk)
- unix shell (tested with bash, dash) and a coreutils implementation (tested with GNU coreutils)


Installation
============

To install project dependencies, run install.sh.


Input data
==========

A certain number of input files are necessary to run the analysis, and must be added manually.

gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx: RepBase repeat counts in hg19 and correlation to A/B profile
gf/Table-AouBouAlways.xlsx: A/B profile for several cell lines
gf/Liste.ProtoSil.Forts-Etudiants-14nov2019.xlsx: Experimental lists of strong protosilencers

Additional CSV files list the source data used in the analysis, which are downloaded automatically.


Usage
=====

To re-run the entire pipeline, including dependency installation and data fetching, run run.sh.


Details
=======

The run.sh script calls the various parts of the pipeline in order and aborts on failure.

- clean.sh:	clean directory from any fetched data, and generated results

- install.sh:	install required packages system-wide

- install.R:	install required R packages and dependencies

- init.R:	make data and result directories and fetch remote data

- prep.R:	convert and filter input data into usable formats

- count.sh:	count single and composite elements along reference genome in 100kb windows

- tabs.R:	generate count summary tables and split data into A/B classes and subclasses

- bedgraph.sh:	generate bedgraph files for count visualization from all count data

- plot.R:	generate basic plots for all count data


Used package versions
=====================

$ bedtools --version
bedtools v2.27.1

> sessionInfo()
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 19.10

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
[1] ggplot2_3.2.1     doParallel_1.0.15 iterators_1.0.12  foreach_1.4.7
[5] readxl_1.3.1      dplyr_0.8.3

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2       magrittr_1.5     tidyselect_0.2.5 munsell_0.5.0
 [5] colorspace_1.4-1 R6_2.4.0         rlang_0.4.0      grid_3.6.1
 [9] gtable_0.3.0     withr_2.1.2      lazyeval_0.2.2   assertthat_0.2.1
[13] tibble_2.1.3     crayon_1.3.4     purrr_0.3.3      codetools_0.2-16
[17] glue_1.3.1       compiler_3.6.1   pillar_1.4.2     cellranger_1.1.0
[21] scales_1.0.0     pkgconfig_2.0.3
