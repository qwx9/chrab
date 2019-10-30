#!/bin/bash -e

getchr(){
	awk '
	NR > 1{
		a[$1]=""
	}
	END{
		s=""
		for(i in a)
			if(s == "")
				s = "\(" i
			else
				s = s "\|" i
		print s "\)"
	}' gf/ab.tsv
}

# tree structure
mkdir -p hg19 huvec imr90 cnt plot

if [[ ! -f gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx || ! -f gf/Table-AouBouAlways.xlsx ]]; then
	echo 'missing input files in gf/ subdirectory!'
	exit 1
fi

# hg19
if [[ ! -e hg19/README ]]; then
	w3m -dump http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ >hg19/README
fi
# hg19 repeatmasker output
wget -nc -P hg19/ 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz'
# hg19 known genes (ncbi refseq)
wget -nc -P hg19/ 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz'
# hg19 ensembl regulation GFF, convert to BED and extract promoters and enhancers only
wget -nc -P hg19/ ftp://ftp.ensembl.org/pub/grch37/release-98/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz
# hg19 %gc track
wget -nc -P hg19/ 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.wigVarStep.gz'
# hg19 chromHMM tracks
wget -nc -P hg19/ 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/*.gz'

# download files specified in csv's, convert excel stuff
Rscript init.R

# hg19 chromosome info (for windowing), needs gf/ab.tsv to restrict chromosome names
if [[ ! -f hg19/hg19.txt ]]; then
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
		-e "select chrom, size from hg19.chromInfo" |\
			grep "^`getchr`\\s" |\
			sort -k1V,1 >hg19/hg19.txt
fi
