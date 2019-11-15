require(dplyr)
require(ggplot2)
require(ggridges)
require(gridExtra)
require(doParallel)

ggviolin <- function(ab, var, class){
	ggplot(ab, aes(!!sym(class), !!sym(var))) +
		geom_violin(na.rm=TRUE) +
		theme(axis.text.x=element_text(angle=90,vjust=-0.1))
}

ggridge <- function(ab, var, class){
	ggplot(ab, aes(!!sym(var), !!sym(class), fill=!!sym(class))) +
		scale_fill_discrete(guide=FALSE) +
		geom_density_ridges2(na.rm=TRUE)
}

ggscatter <- function(ab, var, pc1){
	ggplot(ab, aes(!!sym(var), !!sym(pc1))) +
		geom_hex(na.rm=TRUE)
}

ab <- read.table("tabs/aball.tsv.gz", header=TRUE) %>%
	select(-chr, -start, -end, -AorBvec)
l <- colnames(ab)
l <- l[grep("^(hg19w\\.NA_|hg19w\\.[AB]_|class|HUVEC|IMR90)", l, invert=TRUE)]
f <- paste0("plot/", l, ".pdf")
i <- which(file.access(f, 4) != 0)
if(length(i) == 0)
	return
l <- l[i]
f <- f[i]
abf <- lapply(l, function(x){
	ab %>%
		select(class, HUVEC, !!sym(x)) %>%
		filter(!grepl("^NA_", class))
})
abnf <- lapply(l, function(x){
	ab %>%
		select(classF, HUVECnoflank, !!sym(x)) %>%
		filter(!grepl("^NA_", classF))
})
nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
l <- foreach(i=l, f=f, abf=abf, abnf=abnf, .inorder=FALSE, .packages=c("ggplot2", "ggridges", "gridExtra")) %dopar% {
	g1 <- ggviolin(abf, i, "class")
	g2 <- ggviolin(abnf, i, "classF")
	g3 <- ggridge(abf, i, "class")
	g4 <- ggridge(abnf, i, "classF")
	g5 <- ggscatter(abf, i, "HUVEC")
	g6 <- ggscatter(abnf, i, "HUVECnoflank")
	pdf(f, width=24, height=20)
	print(grid.arrange(g1, g2, g3, g4, g5, g6))
	dev.off()
}
stopCluster(cl)
