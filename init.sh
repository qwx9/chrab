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
# hg19 sequence and repeatmasker output
wget -nc 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa*'
# hg19 promoters (fasta only, not really useful?)
#wget -nc 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/upstream*gz'
# hg19 ensembl regulation GFF, convert to BED and extract promoters and enhancers only
if [[ ! -f hg19.promenh.20180925.bed ]]; then
	wget -O - ftp://ftp.ensembl.org/pub/grch37/release-98/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz \
		| gunzip -c \
		| gff2bed \
		| grep 'promoter:\|enhancer:' \
		> hg19.promenh.20180925.bed
fi
cd ..

# download files specified in csv's, convert excel stuff
Rscript init.R
