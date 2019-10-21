#!/bin/bash
# tree structure
mkdir -p gf doc hg19 huvec imr90 cnt

# documentation
wget -nc https://www.encodeproject.org/documents/a4a6caad-61ab-4d35-820a-409cadce1121/@@download/attachment/RSEM_quantifications_specifications.txt -O doc/rnaseq.tsv.rsem.spec.txt

# hg19
if [[ ! -e hg19/README ]]; then
	w3m -dump http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ >hg19/README
fi
# hg19 repeatmasker output
if [[ ! -f hg19/hg19.fa.out ]]; then
	wget -nc -O hg19/hg19.fa.out.gz 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz'
	gunzip hg19/hg19.fa.out.gz
fi
# hg19 chromosome info (for windowing)
if [[ ! -f hg19/hg19.txt ]]; then
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo" >hg19/hg19.txt
fi
# hg19 ensembl regulation GFF, convert to BED and extract promoters and enhancers only
if [[ ! -f hg19/hg19.promenh.20180925.bed ]]; then
	wget -O - ftp://ftp.ensembl.org/pub/grch37/release-98/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz |\
		gunzip -c |\
		gff2bed |\
		sed -n '/promoter:\|enhancer:/s/^/chr/p' \
		>hg19/hg19.promenh.20180925.bed
fi

# download files specified in csv's, convert excel stuff
Rscript init.R

# extract pro-A repseqs
awk '
FNR == NR && NR > 1{
	if($4 > 0.04)
		s[$3] = ""
}
FNR != NR{
	if($10 in s)
		printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $10)
}
' gf/huvec.repseq.tsv hg19/hg19.fa.out > cnt/huvec.repseq.proa.bed

# extract pro-B repseqs
awk '
FNR == NR && NR > 1{
	if($4 < -0.01)
		s[$3] = ""
}
FNR != NR{
	if($10 in s)
		printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $10)
}
' gf/huvec.repseq.tsv hg19/hg19.fa.out > cnt/huvec.repseq.prob.bed
