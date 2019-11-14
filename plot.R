require(dplyr)
require(ggplot2)
require(doParallel)

ggviolin <- function(ab, var, class){
	f <- paste0("plot/violin.", var, sub("F", ".noflank", substring(class, 6)), ".pdf")
	if(file.access(f, 4) != 0){
		pdf(f, width=12.1, height=9.7)
		g <- ggplot(ab, aes(!!sym(class), !!sym(var))) +
			geom_violin(na.rm=TRUE) +
			theme(axis.text.x=element_text(angle=90,vjust=-0.1))
		print(g)
		dev.off()
	}
}

ggscatter <- function(ab, var, pc1){
	f <- paste0("plot/scatter.", tolower(pc1), ".", var, ".pdf")
	if(file.access(f, 4) != 0){
		pdf(f, width=12.1, height=9.7)
		g <- ggplot(ab, aes(!!sym(var), !!sym(pc1))) +
			geom_hex(na.rm=TRUE)
		print(g)
		dev.off()
	}
}

ab <- read.table("tabs/aball.tsv.gz", header=TRUE) %>%
	select(-chr, -start, -end, -AorBvec)
l <- colnames(ab)
l <- l[grep("^(hg19w\\.NA_|hg19w\\.[AB]_|class|HUVEC|IMR90)", l, invert=TRUE)]
ab <- lapply(l, function(x) select(ab, class, classF, HUVEC, HUVECnoflank, !!sym(x)))
nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
l <- foreach(i=l, ab=ab, .inorder=FALSE, .packages="ggplot2") %dopar% {
	ggviolin(ab, i, "class")
	ggviolin(ab, i, "classF")
	ggscatter(ab, i, "HUVEC")
	ggscatter(ab, i, "HUVECnoflank")
}
stopCluster(cl)
