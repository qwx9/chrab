#!/bin/sh -e
# count elements by 100kb windows along the reference genome
# files in the prep subdirectory will all be counted, as such any files that
# should not be counted (already in 100kb bins or not just intervals) are
# already in the cnt subdirectory

# split hg19 into 100kb windows
bedtools makewindows -g prep/hg19.txt -w 100000 >prep/hg19w.bed

cp huvec/h3k4me3.bed.gz prep/huvec.h3k4me3.bed.gz
cp huvec/h3k27ac.bed.gz prep/huvec.h3k27ac.bed.gz
cp huvec/h3k27me3.bed.gz prep/huvec.h3k27me3.bed.gz
cp huvec/dhs.rep1.bed.gz prep/huvec.dhs.rep1.bed.gz

bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i prep/huvec.groseq.allrep.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	gzip -c >prep/huvec.groseq.mean.bed.gz
bedtools merge -s -c 4,5,6 -o distinct,sum,distinct -i prep/huvec.groseq.allrep.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	gzip -c >prep/huvec.groseq.sum.bed.gz
bedtools map -c 5 -o mean -null 0 -a prep/hg19w.bed -b prep/huvec.groseq.mean.bed.gz |\
	gzip -c >cnt/huvec.groseq.meanofmean.bed.gz
bedtools map -c 5 -o mean -null 0 -a prep/hg19w.bed -b prep/huvec.groseq.sum.bed.gz |\
	gzip -c >cnt/huvec.groseq.meanofsum.bed.gz
bedtools map -c 5 -o sum -null 0 -a prep/hg19w.bed -b prep/huvec.groseq.sum.bed.gz |\
	gzip -c >cnt/huvec.groseq.sumofsum.bed.gz
bedtools merge -s -d 1 -c 4,5,6 -o distinct,distinct,distinct -i prep/huvec.groseq.sum.bed.gz |\
	gzip -c >prep/huvec.groseq.merged.bed.gz

zcat prep/huvec.chromhmm.1.active_promoter.bed.gz prep/huvec.chromhmm.2.weak_promoter.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	bedtools merge -c 4 -o collapse -i - |\
	gzip -c >prep/huvec.chromhmm.any.active.promoters.bed.gz
zcat prep/huvec.chromhmm.1.active_promoter.bed.gz prep/huvec.chromhmm.2.weak_promoter.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	bedtools merge -d 1000 -c 4 -o collapse -i - >/tmp/bedtools.merge
awk '$4 ~ /^[2,]+$/{ print $0 }' /tmp/bedtools.merge |\
	gzip -c >prep/huvec.chromhmm.merged.weak.promoters.bed.gz
awk '$4 !~ /^[2,]+$/{ print $0 }' /tmp/bedtools.merge |\
	gzip -c >prep/huvec.chromhmm.merged.strong.promoters.bed.gz
rm /tmp/bedtools.merge

zcat prep/huvec.chromhmm.4.strong_enhancer.bed.gz prep/huvec.chromhmm.5.strong_enhancer.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	bedtools merge -c 4 -o collapse -i - |\
	gzip -c >prep/huvec.chromhmm.any.strong.enhancers.bed.gz
zcat prep/huvec.chromhmm.6.weak_enhancer.bed.gz prep/huvec.chromhmm.7.weak_enhancer.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	bedtools merge -c 4 -o collapse -i - |\
	gzip -c >prep/huvec.chromhmm.any.weak.enhancers.bed.gz

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd prep
ls -1 *.gz repseq/*.gz |\
	grep -v 'hg19\.gc\|groseq\.\(allrep\.\|mean\.\|sum\.\)' |\
	xargs -I {} -P 8 bash -c \
		'bedtools coverage -counts -a hg19w.bed -b {} | gzip -c > ../cnt/{}'
cp hg19.gc.bed.gz ../cnt/hg19.gc.bed.gz
