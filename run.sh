#!/bin/sh -e
# rerun entire pipeline from scratch
./clean.sh
./install.sh
Rscript init.R
Rscript prep.R
./count.sh
Rscript tabs.R
Rscript lm.R
Rscript plot.R
Rscript bedgraph.R
