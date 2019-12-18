require(dplyr)
require(grid)
require(gridExtra)
require(ggplot2)
require(ggridges)
require(doParallel)

ggviolin <- function(ab, var, class){
	n <- ab %>%
		group_by(!!sym(class)) %>%
		summarize(m=ifelse(typeof(!!sym(var)) == "integer",
			sum(!!sym(var)), rep(NA)))
	n <- setNames(n$m, as.character(pull(n, !!sym(class))))
	ggplot(ab, aes(!!sym(class), !!sym(var))) +
		geom_violin(na.rm=TRUE) +
		theme(axis.text.x=element_text(angle=90, vjust=-0.1)) +
		annotate("text", x=names(n), y=max(ab[,var]) + (max(nchar(n))-1)/.pt,
			label=n, angle=90, size=3)
}

ggridge <- function(ab, var, class){
	ggplot(ab, aes(!!sym(var), !!sym(class), fill=!!sym(class))) +
		scale_fill_discrete(guide=FALSE) +
		geom_density_ridges2(na.rm=TRUE)
}

ggscatter <- function(ab, var, pc1){
	n <- sum(ab[,var])
	Ea <- round(sum(ab[ab$class2=="A",var]) / n * 100, 1)
	EA <- round(sum(ab[ab$class3=="A",var], na.rm=TRUE) / n * 100, 1)
	Eb <- round(sum(ab[ab$class2=="B",var]) / n * 100, 1)
	EB <- round(sum(ab[ab$class3=="B",var], na.rm=TRUE) / n * 100, 1)
	s <- paste0("A Enrichment: ", Ea, "%, AlwaysA: ", EA, "%\n",
		"B Enrichment: ", Eb, "%, AlwaysB: ", EB, "%\n")
	g <- ggplot(ab, aes(!!sym(var), !!sym(pc1))) +
		geom_hex(na.rm=TRUE)
	cnt <- max(layer_data(g, 1)[,"count"])
	if(!is.infinite(cnt)){
		g <- g +
			scale_fill_gradientn(colors=c("black", "blue", "red", "yellow"),
				values=c(0, 1/50, 1/2, 1),
				breaks=unique(round(seq(0, cnt, length.out=12), 0))) +
			guides(fill=guide_colorbar(barheight=30))
	}
	grid.arrange(g, top=textGrob(s))
}

ab <- read.table("tabs/aball.tsv.gz", header=TRUE) %>%
	select(-chr, -start, -end)
l <- colnames(ab)
l2 <- l[grep("repseq\\.only", l)]
l <- l[grep("^(NA_|[AB]_|class|HUVEC|IMR90|huvec\\.repseq\\.only)", l, invert=TRUE)]
f <- c(paste0("plot/", l, ".pdf"), paste0("plot/repseq/", l2, ".pdf"))
l <- c(l, l2)
i <- which(file.access(f, 4) != 0)
if(length(i) == 0)
	quit()
l <- l[i]
f <- f[i]
abf <- ab %>%
	filter(!grepl("^NA_", class)) %>%
	mutate(class=droplevels(class),
		class2=ifelse(grepl("^A_", class), "A", "B"),
		class3=ifelse(grepl("^A_.*AlwaysA", class), "A",
			ifelse(grepl("^B_.*AlwaysB", class), "B",
			NA)))
abf <- lapply(l, function(x){
	abf %>%
		select(class, class2, class3, HUVEC, !!sym(x))
})
abnf <- ab %>%
	filter(!grepl("^NA_", classF)) %>%
	mutate(classF=droplevels(classF),
		class2=ifelse(grepl("^A_", classF), "A", "B"),
		class3=ifelse(grepl("^A_.*AlwaysA", classF), "A",
			ifelse(grepl("^B_.*AlwaysB", classF), "B",
			NA)))
abnf <- lapply(l, function(x){
	abnf %>%
		select(classF, class2, class3, HUVECnoflank, !!sym(x))
})
nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
l <- foreach(i=l, f=f, abf=abf, abnf=abnf, .inorder=FALSE, .multicombine=TRUE, .packages=c("ggplot2", "ggridges", "grid", "gridExtra", "dplyr")) %dopar% {
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
