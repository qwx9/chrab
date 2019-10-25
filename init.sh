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
mkdir -p gf hg19 huvec imr90 cnt plot

# hg19
if [[ ! -e hg19/README ]]; then
	w3m -dump http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/ >hg19/README
fi
# hg19 repeatmasker output
if [[ ! -f hg19/hg19.fa.out ]]; then
	wget -nc -O hg19/hg19.fa.out.gz 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.out.gz'
	gunzip hg19/hg19.fa.out.gz
fi
# hg19 known genes annotation
if [[ ! -f hg19/hg19.refgene.txt.gz ]]; then
	wget -nc -O hg19/hg19.refgene.txt.gz 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
fi
# hg19 ensembl regulation GFF, convert to BED and extract promoters and enhancers only
if [[ ! -f hg19/hg19.promenh.20180925.bed.gz ]]; then
	wget -O - ftp://ftp.ensembl.org/pub/grch37/release-98/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz |\
		gunzip -c |\
		gff2bed |\
		sed -n '/promoter:\|enhancer:/s/^/chr/p' |\
		gzip -c >hg19/hg19.promenh.20180925.bed.gz
fi
# hg19 %gc track
if [[ ! hg19/hg19.gc5Base.wig.gz ]]; then
	wget -nc -O hg19/hg19.gc5base.wig.gz 'ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.gc5Base.wig.gz'
fi

# download files specified in csv's, convert excel stuff
Rscript init.R

# hg19 chromosome info (for windowing), needs gf/ab.tsv to restrict chromosome names
if [[ ! -f hg19/hg19.txt ]]; then
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A \
		-e "select chrom, size from hg19.chromInfo" |\
			grep "^`getchr`\\s" |\
			sort -k1V,1 >hg19/hg19.txt
fi
