require(dplyr)

# generate a/b classes
ab <- read.table("prep/ab.bed", header=TRUE)
ab$ngene <- read.table("cnt/hg19w.hg19.refseq.bed.gz") %>%
	filter(V1 %in% unique(ab$chr)) %>%
	pull(V4)
ab$nact <- read.table("cnt/hg19w.huvec.proa.genes.bed.gz") %>%
	filter(V1 %in% unique(ab$chr)) %>%
	pull(V4)
ab <- ab %>%
	mutate(class1=ifelse(HUVEC < 0, "B", "A"),
		class2=ifelse(ngene >= 4, "highgenedensity", ifelse(ngene > 0, "normalgenedensity", "nogene")),
		class3=ifelse(nact > 0, "hasactive", "noactive"),
		class=paste(class1, class2, AorBvec, class3, sep="_")) %>%
	select(-class1, -class2, -class3)
write.table(ab, "tabs/class.txt", sep="\t", row.names=FALSE, quote=FALSE)
