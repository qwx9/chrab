#!/bin/sh -e

mergenames(){
	gunzip -c $* |\
		awk '
		{
			i=$1$4
			if(e[i]==""){
				g[i]=$4 "\t" $5
				s[i]=$2
				e[i]=$3
				c[i]=$1
				ss[i]=$6
			}else{
				if(ss[i] != $6){
					err="strand conflict!"
					exit -1
				}
				e[i]=$3
				if(s[i] > $2)
					s[i]=$2
			}
		}
		END{
			if(err != ""){
				print err
				exit 2
			}
			for(i in g)
				print c[i] "\t" s[i] "\t" e[i] "\t" g[i] "\t" ss[i]
		}' |\
		sort -k1,1 -k2n,2 |\
		gzip -c
}

cp huvec/h3k4me3.bed.gz prep/huvec.h3k4me3.bed.gz
cp huvec/h3k27ac.bed.gz prep/huvec.h3k27ac.bed.gz
cp huvec/dhs.rep1.bed.gz prep/huvec.dhs.rep1.bed.gz

bedtools merge -s -c 4,5,6 -o distinct,mean,distinct -i prep/huvec.groseq.allrep.bed.gz |\
	gzip -c >prep/huvec.groseq.mean.bed.gz
bedtools merge -s -d 1 -c 4,5,6 -o distinct,distinct,distinct -i prep/huvec.groseq.mean.bed.gz |\
	gzip -c >prep/huvec.groseq.merged.bed.gz
bedtools intersect -s -a prep/hg19.refseq.bed.gz -b prep/huvec.groseq.merged.bed.gz |\
	gzip -c >prep/huvec.refseq.inter.bed.gz
mergenames prep/huvec.refseq.inter.bed.gz >prep/huvec.active.genes.bed.gz

zcat prep/huvec.chromhmm.1.active_promoter.bed.gz prep/huvec.chromhmm.2.weak_promoter.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	bedtools merge -i - |\
	gzip -c >prep/huvec.chromhmm.all.active.promoters.bed.gz
zcat prep/huvec.chromhmm.4.strong_enhancer.bed.gz prep/huvec.chromhmm.5.strong_enhancer.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	bedtools merge -i - |\
	gzip -c >prep/huvec.chromhmm.all.strong.enhancers.bed.gz
zcat prep/huvec.chromhmm.6.weak_enhancer.bed.gz prep/huvec.chromhmm.7.weak_enhancer.bed.gz |\
	sort -k1V,1 -k2n,2 |\
	bedtools merge -i - |\
	gzip -c >prep/huvec.chromhmm.all.weak.enhancers.bed.gz

bedtools intersect -a prep/hg19.prom.gff.gz -b prep/huvec.chromhmm.all.active.promoters.bed.gz |\
	gzip -c >prep/huvec.active.promoters.bed.gz

# split hg19 into 100kb windows
bedtools makewindows -g prep/hg19.txt -w 100000 >prep/hg19w.bed

# count all elements (proa/prob/repseq) per 100kb window along hg19
cd prep
ls -1 *.gz repseq/*.gz |\
	grep -v 'hg19\.gc\|groseq\.allrep' |\
	xargs -I {} -P 8 bash -c \
		'bedtools coverage -counts -a hg19w.bed -b {} | gzip -c > ../cnt/{}'
cp hg19.gc.bed.gz ../cnt/hg19.gc.bed.gz
