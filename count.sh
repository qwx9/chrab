#!/bin/bash -e

bedsub(){
	bedtools subtract -A -a $1 -b ${@:2} |\
	gzip -c
}

bedinter(){
	bedtools intersect -a $1 -b ${@:2} |\
	gzip -c
}

# active promoters
bedinter huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.prm.bed.gz
# active enhancers
bedsub huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.enh.bed.gz
# open promoters/enhancers
bedsub huvec/dhs.rep1.bed.gz huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.open.bed.gz
# inactive promoters/enhancers
bedsub hg19/hg19.promenh.20180925.bed.gz cnt/huvec.proa.prm.bed.gz cnt/huvec.proa.enh.bed.gz cnt/huvec.proa.open.bed.gz >cnt/huvec.prob.prm.enh.bed.gz

bedtools merge -d 1 -i huvec/groseq.rep3.bedGraph.gz |\
	gzip -c >cnt/huvec.groseq.bed.gz
bedinter cnt/huvec.groseq.bed.gz cnt/hg19.refseq.bed.gz >cnt/huvec.proa.genes.bed.gz

# copy individual elements for counting
cp huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz huvec/dhs.rep1.bed.gz hg19/hg19.promenh.20180925.bed.gz cnt/

# split hg19 into 100kb windows
bedtools makewindows -g hg19/hg19.txt -w 100000 >cnt/hg19w.bed

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd cnt
ls -1 *.gz |\
	xargs -I {} -P 8 bash -c 'bedtools coverage -counts -a hg19w.bed -b {} | gzip -c > hg19w.cnt.{}'
