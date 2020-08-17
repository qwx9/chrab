#!/bin/sh -e
# count elements by 100kb windows along the reference genome
# files in the prep subdirectory will all be counted, as such any files that
# should not be counted (already in 100kb bins or not just intervals) are
# already in the cnt subdirectory

# split hg19 into 100kb windows
bedtools makewindows -g prep/hg19.txt -w 100000 >prep/hg19w.bed

for i in huvec gm12878 bcell; do
	for j in h3k4me3 h3k27ac h3k27me3 h3k4me3 dhs.rep1; do
		if test -f $i/$j.bed.gz; then
			cp $i/$j.bed.gz prep/$i.$j.bed.gz
		fi
	done
done
zcat huvec/h3k9me3.bedgraph.gz |\
	grep -v '^\(track\|chrM\)' |\
	gzip -c >prep/huvec.h3k9me3.bed.gz
for i in gm12878 bcell; do
	zcat $i/h3k9me3.bed.gz |\
		sort -k1V,1 -k2n,2 |\
		awk '{print $1 "\t" $2 "\t" $3 "\t" $5}' |\
		gzip -c >prep/$i.h3k9me3.bed.gz
done

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

for i in huvec gm12878 bcell; do
	bedtools map -c 4 -o mean -null 0 -a prep/hg19w.bed -b prep/$i.h3k9me3.bed.gz |\
		gzip -c >cnt/$i.h3k9me3.mean.bed.gz
	bedtools map -c 4 -o sum -null 0 -a prep/hg19w.bed -b prep/$i.h3k9me3.bed.gz |\
		gzip -c >cnt/$i.h3k9me3.sum.bed.gz
	bedtools coverage -a prep/$i.ab.bed -b prep/$i.h3k9me3.bed.gz |\
		awk '{print $1 "\t" $2 "\t" $3 "\t" $10}' |\
		gzip -c >cnt/$i.h3k9me3.coverage.bed.gz

	zcat prep/$i.chromhmm.1.active_promoter.bed.gz prep/$i.chromhmm.2.weak_promoter.bed.gz |\
		sort -k1V,1 -k2n,2 |\
		bedtools merge -c 4 -o collapse -i - |\
		gzip -c >prep/$i.chromhmm.any.active.promoters.bed.gz
	zcat prep/$i.chromhmm.1.active_promoter.bed.gz prep/$i.chromhmm.2.weak_promoter.bed.gz |\
		sort -k1V,1 -k2n,2 |\
		bedtools merge -d 1000 -c 4 -o collapse -i - >/tmp/bedtools.merge
	awk '$4 ~ /^[2,]+$/{ print $0 }' /tmp/bedtools.merge |\
		gzip -c >prep/$i.chromhmm.merged.weak.promoters.bed.gz
	awk '$4 !~ /^[2,]+$/{ print $0 }' /tmp/bedtools.merge |\
		gzip -c >prep/$i.chromhmm.merged.strong.promoters.bed.gz
	rm /tmp/bedtools.merge

	zcat prep/$i.chromhmm.4.strong_enhancer.bed.gz prep/$i.chromhmm.5.strong_enhancer.bed.gz |\
		sort -k1V,1 -k2n,2 |\
		bedtools merge -c 4 -o collapse -i - |\
		gzip -c >prep/$i.chromhmm.any.strong.enhancers.bed.gz
	if test -f prep/$i.chromhmm.7.weak_enhancer.bed.gz; then
		zcat prep/$i.chromhmm.6.weak_enhancer.bed.gz prep/$i.chromhmm.7.weak_enhancer.bed.gz |\
			sort -k1V,1 -k2n,2 |\
			bedtools merge -c 4 -o collapse -i - |\
			gzip -c >prep/$i.chromhmm.any.weak.enhancers.bed.gz
	fi
done

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd prep
ls -1 *.gz repseq/*.gz |\
	grep -v 'hg19\.gc\|h3k9me3\|groseq\.\(allrep\.\|mean\.\|sum\.\)' |\
	xargs -I {} -P 8 bash -c \
		'bedtools coverage -counts -a hg19w.bed -b {} | gzip -c > ../cnt/{}'
cp hg19.gc.bed.gz ../cnt/hg19.gc.bed.gz
