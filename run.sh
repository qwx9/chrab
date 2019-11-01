#!/bin/sh -e
./clean.sh
./install.sh
./installtools.sh
Rscript init.R
Rscript prep.R
./count.sh
Rscript tabs.R
./bedgraph.sh
Rscript plot.R
