Dependencies
============

Development and testing was done on an up-to-date Ubuntu Linux 19.10 instance.

Requirements:
- bash: shell scripts
- R >= 3.6.0: R scripts
- bedops: gff to bed conversion
- w3m: html to plain text conversion


Installation
============

To install project dependencies, run install.sh.
To install development tools, run installtools.sh.


Input data
==========

A certain number of input files are necessary to run the analysis.

encode.csv: files to be exported from ENCODE, later used for counting (provided)
gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx: RepBase repeat counts in hg19 and correlation to A/B profile
gf/Table-AouBouAlways.xlsx: A/B profile for several cell lines


Usage
=====

To re-run the entire pipeline, including dependency installation and data fetching, run run.sh.

To only create the directory tree and download/convert used data, run init.sh.

To only count elements along the reference genome, run count.sh.

To only create the plots based on the count dta, run plot.R using Rscript.


Files
=====

*.R: R scripts for importing data, returning results and figures
*.sh: bash scripts


Used package versions
=====================

$ bedtools --version
bedtools v2.27.1

$ bedops --version | grep version
  version:  2.4.36 (typical)

$ R --version | sed 1q
R version 3.6.1 (2019-07-05) -- "Action of the Toes"

FIXME: R sessionInfo() (packages)
