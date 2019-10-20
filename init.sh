#!/bin/sh
# tree structure
mkdir -p doc hg19 huvec imr90

# documentation
wget https://www.encodeproject.org/documents/a4a6caad-61ab-4d35-820a-409cadce1121/@@download/attachment/RSEM_quantifications_specifications.txt -O doc/rnaseq.tsv.rsem.spec.txt

# genome reference files
cd hg19
w3m -dump http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ >README
wget --timestamping 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa*'
