#!/bin/bash
# active promoters
bedtools intersect -a huvec/h3k4me3.bed.gz -b huvec/h3k27ac.bed.gz >cnt/huvec.proa.prm.bed
# active enhancers
bedtools subtract -A -a huvec/h3k4me3.bed.gz -b huvec/h3k27ac.bed.gz >cnt/huvec.proa.enh.bed
# open promoters/enhancers
bedtools intersect -a huvec/dhs.rep1.bed.gz -b huvec/dhs.rep2.bed.gz >cnt/huvec.dhs.bed
bedtools subtract -A -a cnt/huvec.dhs.bed -b huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.proa.open.bed
# pro-A repseqs
awk '
FNR == NR && NR > 1{
	if($4 > 0.04)
		s[$3] = ""
}
FNR != NR{
	if($10 in s)
		printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $10)
}
' gf/huvec.repseq.tsv hg19/hg19.fa.out >cnt/huvec.proa.repseq.bed

# huvec pro-B
bedtools subtract -A -a hg19/hg19.promenh.20180925.bed -b cnt/huvec.proa.prm.bed cnt/huvec.proa.open.bed >cnt/huvec.prob.prm.enh.bed
# pro-B repseqs
awk '
FNR == NR && NR > 1{
	if($4 < -0.01)
		s[$3] = ""
}
FNR != NR{
	if($10 in s)
		printf("%s\t%s\t%s\t%s\n", $5, $6, $7, $10)
}
' gf/huvec.repseq.tsv hg19/hg19.fa.out >cnt/huvec.prob.repseq.bed

# split hg19 into 100kb windows
bedtools makewindows -g hg19/hg19.txt -w 100000 >cnt/hg19w.txt

# count all elements per 100kb window along hg19
cd cnt
for i in huvec.pro{a,b}.*.bed; do
	bedtools coverage -counts -a hg19w.txt -b $i >hg19w.cnt.$i
done
