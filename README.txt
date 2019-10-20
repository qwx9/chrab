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

To create the directory tree and download/convert used data, run init.sh.


Files
=====

encode.csv: files exported from ENCODE
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
