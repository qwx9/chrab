#!/bin/bash
# huvec pro-A
bedtools intersect -a huvec/h3k4me3.bed.gz -b huvec/h3k27ac.bed.gz >cnt/huvec.prom.actif.bed
bedtools subtract -A -a huvec/h3k4me3.bed.gz -b huvec/h3k27ac.bed.gz >cnt/huvec.enh.actif.bed
bedtools intersect -a huvec/dhs.rep1.bed.gz -b huvec/dhs.rep2.bed.gz >cnt/huvec.dhs.bed
bedtools subtract -A -a cnt/huvec.dhs.bed -b huvec/h3k4me3.bed.gz huvec/h3k27ac.bed.gz >cnt/huvec.prom.enh.open.bed

# huvec pro-B
bedtools subtract -A -a hg19/hg19.promenh.20180925.bed -b cnt/huvec.prom.actif.bed cnt/huvec.prom.enh.open.bed >cnt/huvec.prom.enh.inactif.bed
