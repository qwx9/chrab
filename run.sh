#!/bin/bash
./clean.sh
./install.sh
./installtools.sh
Rscript init.R
Rscript prep.R
./count.sh
Rscript plot.R
