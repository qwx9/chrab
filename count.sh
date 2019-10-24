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
	' gf/huvec.repseq.tsv hg19/hg19.fa.out
}

getrep(){
	f=huvec.repseq.only.`echo $1 | sed 's,/,-,g;s,(\|),,g'`.bed
	grep "$1\s" hg19/hg19.fa.out |\
		awk '{printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $10)}' >cnt/$f
}

# xargs black magic
reps2beds(){
	export -f getrep
	awk 'NR>1{print $3}' gf/huvec.repseq.tsv |\
		xargs -P 8 -I {} bash -c 'getrep "$@"' _ {}
}

# hg19 known genes
gunzip -c hg19/hg19.refgene.txt.gz |\
	awk '{printf("%s\t%s\t%s\t%s\n", $3, $5, $6, $4)}' |\
	sort -k 1d,1 |\
	uniq -d >cnt/hg19.refgene.bed

# active genes list
gunzip -c huvec/groseq.txt.gz |\
	awk 'NR>1{printf("%s\t%s\t%s\t%s\n", $2, $3, $4, $5)}' |\
	sort -k 1d,1 |\
	uniq -d >cnt/huvec.proa.genes.bed

# repseqs
rep2bed '$1 == "LTR" || $1 == "LINE"' >cnt/huvec.repseq.ltrline.bed
rep2bed '$1 == "LTR"' >cnt/huvec.repseq.ltr.bed
rep2bed '$1 == "LINE" && $2 == "L1"' >cnt/huvec.repseq.l1.bed
rep2bed '$4 < -0.01' >cnt/huvec.prob.repseq.bed
rep2bed '$4 < -0.01 && ($1 == "LTR" || $1 == "LINE")' >cnt/huvec.prob.repseq.ltrline.bed
rep2bed '$4 < -0.01 && $1 == "LTR"' >cnt/huvec.prob.repseq.ltr.bed
rep2bed '$4 < -0.01 && $1 == "LINE"' >cnt/huvec.prob.repseq.line.bed
rep2bed '$4 < -0.01 && ($1 == "Simple_repeat" || $1 == "Low_complexity" || $1 == "Satellite")' >cnt/huvec.prob.repseq.sr.lc.sat.bed
rep2bed '$4 > 0.04' >cnt/huvec.proa.repseq.bed

# individual repseqs
reps2beds

# active promoters
bedtools intersect -a huvec/h3k4me3.bed.gz -b huvec/h3k27ac.bed.gz >cnt/huvec.proa.prm.bed
# active enhancers
bedtools subtract -A -a huvec/h3k4me3.bed.gz -b huvec/h3k27ac.bed.gz >cnt/huvec.proa.enh.bed
# open promoters/enhancers
bedtools subtract -A -a huvec/dhs.rep1.bed.gz -b huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.open.bed
# inactive promoters/enhancers
bedtools subtract -A -a hg19/hg19.promenh.20180925.bed -b cnt/huvec.proa.prm.bed cnt/huvec.proa.enh.bed cnt/huvec.proa.open.bed >cnt/huvec.prob.prm.enh.bed

# split hg19 into 100kb windows
bedtools makewindows -g hg19/hg19.txt -w 100000 >cnt/hg19w.txt

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd cnt
ls -1 hg19.refgene.bed huvec.pro{a,b}.*.bed huvec.repseq.*.bed |\
	xargs -I {} -P 8 bash -c 'bedtools coverage -counts -a hg19w.txt -b {} > hg19w.cnt.{}'
