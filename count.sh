#!/bin/bash -e

bedsub(){
	bedtools subtract -A -a $1 -b ${@:2} |\
	gzip -c
}

bedinter(){
	bedtools intersect -a $1 -b ${@:2} |\
	gzip -c
}

cp huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz huvec/dhs.rep1.bed.gz prep/
# active promoters
bedinter prep/h3k4me3.bed.gz prep/h3k27ac.bed.gz >prep/huvec.proa.prm.bed.gz
# active enhancers
bedsub prep/h3k4me3.bed.gz prep/h3k27ac.bed.gz >prep/huvec.proa.enh.bed.gz
# open promoters/enhancers
bedsub prep/dhs.rep1.bed.gz prep/h3k4me3.bed.gz prep/h3k27ac.bed.gz >prep/huvec.proa.open.bed.gz
# inactive promoters/enhancers
bedsub prep/hg19.promenh.gff.gz prep/huvec.proa.prm.bed.gz prep/huvec.proa.enh.bed.gz prep/huvec.proa.open.bed.gz >prep/huvec.prob.prm.enh.bed.gz

bedtools merge -d 1 -i huvec/groseq.rep3.bedGraph.gz |\
	gzip -c >prep/huvec.groseq.bed.gz
bedinter prep/huvec.groseq.bed.gz prep/hg19.refseq.bed.gz >prep/huvec.proa.genes.bed.gz

# split hg19 into 100kb windows
bedtools makewindows -g prep/hg19.txt -w 100000 >prep/hg19w.bed

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd prep
ls -1 *.gz |\
	xargs -I {} -P 8 bash -c \
		'bedtools coverage -counts -a hg19w.bed -b {} | gzip -c > ../cnt/hg19w.{}'
