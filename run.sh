#!/bin/bash
./clean.sh
./install.sh
./installtools.sh
./init.sh
./count.sh
Rscript plot.R
