#!/bin/sh -e
# rerun entire pipeline from scratch
exec >/dev/null
./clean.sh
./install.sh
echo '-> downloading input files' >&2
Rscript init.R
echo '-> preparing input files' >&2
Rscript prep.R
echo '-> counting elements' >&2
./count.sh
echo '-> generating count tables' >&2
Rscript tabs.R
echo '-> generating models' >&2
Rscript lm.R
echo '-> generating result plots' >&2
Rscript plot.R
echo '-> generating count bedgraphs' >&2
Rscript bedgraph.R
