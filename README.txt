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

- bedgraph.R:	generate bedgraph files for count visualization from all count data

- plot.R:	generate basic plots for all count data

- score.R:	generate predictive models and bedgraph and plots of predictions


Additional documentation
========================

The matmeth.md file is a detailed overview of the methods applied and choices made in this pipeline.


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
[1] grid      parallel  stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] tidyr_1.0.0       readxl_1.3.1      lmtest_0.9-37     zoo_1.8-6
 [5] gridExtra_2.3     ggridges_0.5.1    ggplot2_3.2.1     dplyr_0.8.3
 [9] doParallel_1.0.15 iterators_1.0.12  foreach_1.4.7

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.3       magrittr_1.5     tidyselect_0.2.5 munsell_0.5.0
 [5] lattice_0.20-38  colorspace_1.4-1 R6_2.4.1         rlang_0.4.2
 [9] plyr_1.8.5       gtable_0.3.0     withr_2.1.2      lazyeval_0.2.2
[13] assertthat_0.2.1 tibble_2.1.3     lifecycle_0.1.0  crayon_1.3.4
[17] purrr_0.3.3      vctrs_0.2.0      codetools_0.2-16 zeallot_0.1.0
[21] glue_1.3.1       cellranger_1.1.0 compiler_3.6.1   pillar_1.4.2
[25] backports_1.1.5  scales_1.1.0     pkgconfig_2.0.3


License
=======

The code in this repository is covered under the MIT license, reviewable in LICENSE.txt


Contributors
============

This project is the result of work by the following people:

Project leader:
Geneviève Fourel <genevieve.fourel@ens-lyon.fr>, Research director at INSERM, France

Programmers:
Konstantinn Bonnet <konstantinn.bonnet@etu.univ-lyon1.fr>, 1st year student in the Bioinformatics Master's Degree of Lyon 1 university, France
Théophile Boyer <theophile.boyer@etu.univ-lyon1.fr>, 1st year student in the Bioinformatics Master's Degree of Lyon 1 university, France

Preliminary work:
Raphael Mourad <raphael.mourad@ibcg.biotoul.fr>, Assistant professor at University of Toulouse III, France

Additional help:
Jean-Baptiste Claude <jean-baptise.claude@ens-lyon.fr>, Bioinformatics Research Engineer, LBMC, Ens Lyon, France
