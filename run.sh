#!/bin/bash -e
./clean.sh
./install.sh
./installtools.sh
Rscript init.R
Rscript prep.R
./count.sh
Rscript tabs.R
Rscript plot.R
