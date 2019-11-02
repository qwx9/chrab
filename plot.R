require(dplyr)
require(ggplot2)

ggviolin <- function(ab, var){
	f <- paste0("plot/violin.", var, ".pdf")
	if(file.access(f, 4) == 0)
		return
	pdf(f, width=12.1, height=9.7)
	g <- ggplot(ab, aes(class, !!sym(var))) +
		geom_violin(na.rm=TRUE) +
		theme(axis.text.x=element_text(angle=90,vjust=-0.1))
	print(g)
	dev.off()
}

ab <- read.table("tabs/aball.tsv.gz", header=TRUE)
ggviolin(ab, "ngene")
ggviolin(ab, "nact")

#ggplot(ab, aes(ngene, HUVEC)) + geom_hex()
#ggplot(ab, aes(nact, HUVEC)) + geom_hex()
#ab$prma <- read.table("cnt/hg19w.cnt.huvec.proa.prm.bed")$V4
#ggplot(ab, aes(prma, HUVEC)) + geom_hex()
