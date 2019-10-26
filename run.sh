#!/bin/bash -e
./install.sh
./installtools.sh
./init.sh
./count.sh
Rscript plot.R
