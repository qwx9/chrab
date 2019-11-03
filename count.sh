#!/bin/sh -e

bedsub(){
	f=$1
	shift
	bedtools subtract -A -a $f -b $* |\
		gzip -c
}

bedinter(){
	f=$1
	shift
	bedtools intersect -a $f -b $* |\
		gzip -c
}

cp huvec/h3k4me3.bed.gz prep/huvec.h3k4me3.bed.gz
cp huvec/h3k27ac.bed.gz prep/huvec.h3k27ac.bed.gz
cp huvec/dhs.rep1.bed.gz prep/huvec.dhs.rep1.bed.gz

# active promoters
bedinter prep/huvec.h3k4me3.bed.gz prep/huvec.h3k27ac.bed.gz >prep/huvec.proa.prm.bed.gz
# active enhancers
bedsub prep/huvec.h3k4me3.bed.gz prep/huvec.h3k27ac.bed.gz >prep/huvec.proa.enh.bed.gz
# open promoters/enhancers
bedsub prep/huvec.dhs.rep1.bed.gz prep/huvec.h3k4me3.bed.gz prep/huvec.h3k27ac.bed.gz >prep/huvec.proa.open.bed.gz
# inactive promoters/enhancers
bedsub prep/hg19.promenh.gff.gz prep/huvec.proa.prm.bed.gz prep/huvec.proa.enh.bed.gz prep/huvec.proa.open.bed.gz >prep/huvec.prob.prm.enh.bed.gz

bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i prep/huvec.groseq.allrep.bed.gz |\
	gzip -c >prep/huvec.groseq.mean.bed.gz
bedtools merge -s -d 1 -c 6 -o distinct -i prep/huvec.groseq.mean.bed.gz |\
	gzip -c >prep/huvec.groseq.merged.bed.gz
bedinter prep/huvec.groseq.merged.bed.gz prep/hg19.refseq.bed.gz >prep/huvec.proa.genes.bed.gz

# split hg19 into 100kb windows
bedtools makewindows -g prep/hg19.txt -w 100000 >prep/hg19w.bed

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd prep
ls -1 *.gz |\
	grep -v 'hg19w\.gc5base\|groseq\.allrep' |\
	xargs -I {} -P 8 bash -c \
		'bedtools coverage -counts -a hg19w.bed -b {} | gzip -c > ../cnt/hg19w.{}'
cp prep/hg19w.hg19.gc5base.bed.gz cnt/
