require(dplyr)
require(ggplot2)

ab <- read.table("tabs/class.txt", header=TRUE)

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
