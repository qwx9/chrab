#!/bin/bash

rep2bed(){
	awk '
	FNR == NR && NR > 1{
		if('"$1"')
			s[$3] = ""
	}
	FNR != NR{
		if($10 in s)
			printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $10)
	}
	' gf/huvec.repseq.tsv hg19/hg19.fa.out |\
	gzip -c
}

getrep(){
	f=huvec.repseq.only.`echo $1 | sed 's,/,-,g;s,(\|),,g'`.bed.gz
	grep "$1\s" hg19/hg19.fa.out |\
		awk '{printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $10)}' |\
		gzip -c >cnt/$f
}

# xargs black magic
reps2beds(){
	export -f getrep
	awk 'NR>1{print $3}' gf/huvec.repseq.tsv |\
		xargs -P 8 -I {} bash -c 'getrep "$@"' _ {}
}

bedsub(){
	bedtools subtract -A -a $1 -b ${@:2} |\
	gzip -c
}

bedinter(){
	bedtools intersect -a $1 -b ${@:2} |\
	gzip -c
}

# hg19 known genes
gunzip -c hg19/hg19.refgene.txt.gz |\
	awk '{printf("%s\t%s\t%s\t%s\n", $3, $5, $6, $4)}' |\
	sort -k 1d,1 |\
	uniq -d |\
	gzip -c >cnt/hg19.refgene.bed.gz

# active genes list
gunzip -c huvec/groseq.txt.gz |\
	awk 'NR>1{printf("%s\t%s\t%s\t%s\n", $2, $3, $4, $5)}' |\
	sort -k 1d,1 |\
	uniq -d |\
	gzip -c >cnt/huvec.proa.genes.bed.gz

# repseqs
rep2bed '$1 == "LTR" || $1 == "LINE"' >cnt/huvec.repseq.ltrline.bed.gz
rep2bed '$1 == "LTR"' >cnt/huvec.repseq.ltr.bed.gz
rep2bed '$1 == "LINE" && $2 == "L1"' >cnt/huvec.repseq.l1.bed.gz
rep2bed '$4 < -0.01' >cnt/huvec.prob.repseq.bed.gz
rep2bed '$4 < -0.01 && ($1 == "LTR" || $1 == "LINE")' >cnt/huvec.prob.repseq.ltrline.bed.gz
rep2bed '$4 < -0.01 && $1 == "LTR"' >cnt/huvec.prob.repseq.ltr.bed.gz
rep2bed '$4 < -0.01 && $1 == "LINE"' >cnt/huvec.prob.repseq.line.bed.gz
rep2bed '$4 < -0.01 && ($1 == "Simple_repeat" || $1 == "Low_complexity" || $1 == "Satellite")' >cnt/huvec.prob.repseq.sr.lc.sat.bed.gz
rep2bed '$4 > 0.04' >cnt/huvec.proa.repseq.bed.gz

# individual repseqs
reps2beds

# active promoters
bedinter huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.prm.bed.gz
# active enhancers
bedsub huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.enh.bed.gz
# open promoters/enhancers
bedsub huvec/dhs.rep1.bed.gz huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.open.bed.gz
# inactive promoters/enhancers
bedsub hg19/hg19.promenh.20180925.bed.gz cnt/huvec.proa.prm.bed.gz cnt/huvec.proa.enh.bed.gz cnt/huvec.proa.open.bed.gz >cnt/huvec.prob.prm.enh.bed.gz

# copy individual elements for counting
cp huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz huvec/dhs.rep1.bed.gz hg19/hg19.promenh.20180925.bed.gz cnt/

# split hg19 into 100kb windows
bedtools makewindows -g hg19/hg19.txt -w 100000 >cnt/hg19w.bed

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd cnt
ls -1 *.gz |\
	xargs -I {} -P 8 bash -c 'bedtools coverage -counts -a hg19w.bed -b {} | gzip -c > hg19w.cnt.{}.gz'
