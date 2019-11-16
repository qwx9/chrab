#!/bin/sh -e
./clean.sh
./install.sh
Rscript init.R
Rscript prep.R
./count.sh
Rscript tabs.R
Rscript bedgraph.R
Rscript plot.R
