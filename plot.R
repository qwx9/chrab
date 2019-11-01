require(dplyr)
require(ggplot2)

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
write.table(ab, "plot/class.txt", row.names=FALSE, quote=FALSE, sep="\t")

# per class violin plots
pdf("plot/violin.ngene.pdf", width=12.1, height=9.7)
	g <- ggplot(ab, aes(class, ngene)) +
		geom_violin(na.rm=TRUE) +
		theme(axis.text.x=element_text(angle=90,vjust=-0.1))
	print(g)
dev.off()
pdf("plot/violin.nact.pdf", width=12.1, height=9.7)
	g <- ggplot(ab, aes(class, ngene)) +
		geom_violin(na.rm=TRUE) +
		theme(axis.text.x=element_text(angle=90,vjust=-0.1))
	print(g)
dev.off()

#ggplot(ab, aes(ngene, HUVEC)) + geom_hex()
#ggplot(ab, aes(nact, HUVEC)) + geom_hex()
#ab$prma <- read.table("cnt/hg19w.cnt.huvec.proa.prm.bed")$V4
#ggplot(ab, aes(prma, HUVEC)) + geom_hex()
