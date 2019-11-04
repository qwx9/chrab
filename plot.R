require(dplyr)
require(ggplot2)

ggviolin <- function(ab, var){
	f <- paste0("plot/violin.", var, ".pdf")
	if(file.access(f, 4) != 0){
		pdf(f, width=12.1, height=9.7)
		g <- ggplot(ab, aes(class, !!sym(var))) +
			geom_violin(na.rm=TRUE) +
			theme(axis.text.x=element_text(angle=90,vjust=-0.1))
		print(g)
		dev.off()
	}
}

ggscatter <- function(ab, var, pc1){
	f <- paste0("plot/scatter.", var, ".pdf")
	if(file.access(f, 4) != 0){
		pdf(f, width=12.1, height=9.7)
		g <- ggplot(ab, aes(!!sym(var), !!sym(pc1))) +
			geom_hex(na.rm=TRUE)
		print(g)
		dev.off()
	}
}

ab <- read.table("tabs/aball.tsv.gz", header=TRUE)
ggviolin(ab, "ngene")
ggviolin(ab, "nact")
for(i in colnames(ab)[-c(1:6, 9)])
	ggscatter(ab, i, "HUVEC")
