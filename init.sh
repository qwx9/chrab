#!/bin/bash
# tree structure
mkdir -p gf doc hg19 huvec imr90

# documentation
wget -nc https://www.encodeproject.org/documents/a4a6caad-61ab-4d35-820a-409cadce1121/@@download/attachment/RSEM_quantifications_specifications.txt -O doc/rnaseq.tsv.rsem.spec.txt

# hg19
if [[ ! -e hg19/README ]]; then
	w3m -dump http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ >hg19/README
fi
cd hg19
wget -nc 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa*'
cd ..

# encode files
Rscript init.R
