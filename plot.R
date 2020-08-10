# generate plots for all counts
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggridges)
library(doParallel)
})
source("lib.R")

# quantile-quantile plot against normal distribution
ggqnorm <- function(d){
	ggplot(d, aes(sample=res)) +
		stat_qq_line(color="blue") +
		stat_qq(size=0.7) +
		ggtitle("Quantile-quantile plot of residuals") +
		xlab("Theoretical normal distribution quantiles") +
		ylab("Standardized residuals")
}

ggfit <- function(d, var, title, method="auto"){
	ggplot(d, aes(fitted, !!sym(var))) +
		geom_point(na.rm=TRUE, color="darkblue", alpha=0.2, shape=20) +
		geom_smooth(method=method, na.rm=TRUE, se=FALSE, color="red", size=0.4) +
		ggtitle(title)
}

# generate correlation heatmap of params in provided data.frame
ggcor <- function(cm){
	# remove lower matrix triangle since it's redundant
	cm[lower.tri(cm)] <- NA
	# convert matrix to data.frame, then use gather to make a tidy
	# data.frame suitable for ggplot2, and return correlation heatmap
	# ggplot object
	cm %>%
		as.data.frame %>%
		mutate(grp=factor(row.names(.), levels=row.names(.))) %>%
		gather(key, v, -grp, na.rm=TRUE, factor_key=TRUE) %>%
		ggplot(aes(key, grp, fill=v, label=round(v,2))) +
			geom_tile(color="white") +
			scale_fill_gradient2(low="blue", high="red", mid="white",
				midpoint=0, limit=c(-1,1), space="Lab",
				name="Pearson\nCorrelation") +
			theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1)) +
				coord_fixed() +
				geom_text(color="black", size=2)
}

# violin plots: takes a data.frame with class and parameter; parameter name;
# class variable to use
ggviolin <- function(ab, var, class){
	m <- ab %>%
		# for non-integral "counts" like %gc, don't label violin plots
		# otherwise, annotate each violin with its number of observations
		group_by(!!sym(class)) %>%
		summarize(m=ifelse(typeof(!!sym(var)) == "integer",
			sum(!!sym(var)), rep(NA)))
	ab %>%
		ggplot(aes(!!sym(class), !!sym(var))) +
			geom_violin(scale="count", trim=TRUE,
				size=0.3, na.rm=TRUE) +
			geom_jitter(size=0.1, height=0, alpha=0.05, show.legend=FALSE) +
			scale_fill_viridis_d(guide=FALSE) +
			theme(axis.text.x=element_text(angle=90, vjust=-0.1)) +
			geom_text(data=m, aes(!!sym(class),
				max(ab[,var], na.rm=TRUE) + (max(nchar(m))-1)/.pt,
				label=m),
				color="black", angle=45, size=3)
}

# ridgeline plots: takes a data.frame with class and parameter; parameter name;
# class variable to use
ggridge <- function(ab, var, class){
	ab %>%
		mutate_at(class, ~factor(!!sym(class), levels=rev(levels(!!sym(class))), ordered=TRUE)) %>%
		ggplot(aes(!!sym(var), !!sym(class), fill=!!sym(class))) +
			geom_density_ridges2(jittered_points=TRUE, size=0.3, alpha=0.6,
				point_size=0.3, point_color="black", point_alpha=0.1,
				rel_min_height=0.005, na.rm=TRUE) +
			scale_fill_ordinal(guide=FALSE)
}

# scatter plots: takes a data.frame with eigenvector and parameter; parameter
# name; eigenvector variable to use
ggscatter <- function(ab, var, pc1){
	# total count across everything
	n <- sum(ab[,var])
	# enrichment factor: proportion in anyA, alwaysA, anyB, alwaysB
	Ea <- round(sum(ab[ab$class2=="A",var]) / n * 100, 1)
	EA <- round(sum(ab[ab$class3=="A",var], na.rm=TRUE) / n * 100, 1)
	Eb <- round(sum(ab[ab$class2=="B",var]) / n * 100, 1)
	EB <- round(sum(ab[ab$class3=="B",var], na.rm=TRUE) / n * 100, 1)
	# construct annotation string
	s <- paste0("A Enrichment: ", Ea, "%, AlwaysA: ", EA, "%\n",
		"B Enrichment: ", Eb, "%, AlwaysB: ", EB, "%\n")
	# initial plot
	g <- ggplot(ab, aes(!!sym(var), !!sym(pc1))) +
		geom_hex(na.rm=TRUE)
	# extract count maximum for any pc1/frequency combination
	cnt <- max(layer_data(g, 1)[,"count"])
	# if there are no counts (e.g. they were filtered out), cnt is -Inf,
	# just paste an empty plot
	# otherwise, add a custom gradient with a special color for very low
	# values and add more ticks for the gradient label
	if(!is.infinite(cnt)){
		g <- g +
			scale_fill_gradientn(colors=c("#cccccc", "blue", "red", "yellow"),
				values=c(0, 1/50, 1/2, 1),
				breaks=unique(round(seq(0, cnt, length.out=12), 0))) +
			guides(fill=guide_colorbar(barheight=30))
	}
	# put enrichment summary on top of scatter plot
	grid.arrange(g, top=textGrob(s))
}

# stat_density_2d uses MASS::kde2d which in turn uses MASS:bandwidth.nrd,
# which subtracts 3rd quantile from 1st for a bandwidth estimate.
# this obviously fails if we have too many 0 counts.
ggdensity2d <- function(ab, var, pc1){
	if(quantile(ab[,var], 0.5) == 0)
		return(NULL)
	ggplot(ab, aes(!!sym(var), !!sym(pc1), color=classc, fill=classc)) +
		stat_density_2d(geom="polygon", alpha=0.05, na.rm=TRUE) +
		scale_color_brewer(name="subclasses", palette="Paired", aesthetics=c("fill", "color"))
}

ggdensity2dxor <- function(ab, var, pc1, name){
	if(quantile(ab[,var], 0.5) == 0)
		return(NULL)
	ab %>%
		mutate(classc=ifelse(classc==name, name, "others")) %>%
		ggplot(aes(!!sym(var), !!sym(pc1), color=classc, fill=classc)) +
			stat_density_2d(geom="polygon", alpha=0.05, na.rm=TRUE) +
			scale_color_manual(name="subclasses", values=c("red", "grey"), aesthetics=c("fill", "color"))
}

# split bins into classes for exploded scatterplots as specified
scattersplit <- function(x){
	cl <- c("1_AlwaysA.HG",
		"2_AorB.HG",
		"3_AlwaysA.!HG.Ac",
		"4_AlwaysA.!HG.!Ac",
		"5_AorB.NG.Ac",
		"6_AorB.NG.!Ac",
		"7_AorB.0G.Ac",
		"8_AorB.0G.!Ac",
		"9_AlwaysB.NG",
		"10_AlwaysB.0G",
		"11_AnyB.HG"
	)
	factor(ifelse(grepl("highgene.*AlwaysA", x), cl[1],
		ifelse(grepl("^A.*highgene.*AorB", x), cl[2],
		ifelse(grepl("[ol]gene.*AlwaysA.*hasact", x), cl[3],
		ifelse(grepl("[ol]gene.*AlwaysA.*noact", x), cl[4],
		ifelse(grepl("lgene.*AorB.*hasact", x), cl[5],
		ifelse(grepl("lgene.*AorB.*noact", x), cl[6],
		ifelse(grepl("nogene.*AorB.*hasact", x), cl[7],
		ifelse(grepl("nogene.*AorB.*noact", x), cl[8],
		ifelse(grepl("lgene.*AlwaysB", x), cl[9],
		ifelse(grepl("nogene.*AlwaysB", x), cl[10],
		ifelse(grepl("^B.*highgene", x), cl[11],
		NA))))))))))),
		levels=cl, ordered=TRUE)
}

# read table of all counts, removing useless columns
ab <- read.table("tabs/aball.tsv.gz", header=TRUE) %>%
	select(-chr, -start, -end)

# set global theme options
th <- theme_classic() +
	theme(panel.grid.major=element_line(size=0.5, linetype="solid", colour="#eeeeee"),
	panel.grid.minor=element_line(size=0.25, linetype="solid", colour="#eeeeee"))
theme_set(th)

# linear models plots
# get epigenomic and genomic parameter lists
l <- lapply(c("lm.eparm.tsv", "lm.gparm.tsv"), read.parms)
# generate correlation heatmap for all parameters
if(file.access("lm/cor.pdf", 4) != 0){
	read.table("tabs/lmparmscor.tsv", header=TRUE) %>%
		ggcor %>%
		write.pdf("lm/cor.pdf")
}

# reorder classes for plots
l <- levels(ab$class)
l <- c(l[grep("^A_nogene", l)], l[grep("^A_high", l)], l[grep("^A_normal", l)],
	l[grep("^B_normal", l)], l[grep("^B_nogene", l)], l[grep("^B_high", l)],
	l[grep("^NA", l)])
ab$class <- factor(ab$class, levels=l, ordered=TRUE)
l <- levels(ab$classF)
l <- c(l[grep("^A_nogene", l)], l[grep("^A_high", l)], l[grep("^A_normal", l)],
	l[grep("^B_normal", l)], l[grep("^B_nogene", l)], l[grep("^B_high", l)],
	l[grep("^NA", l)])
ab$classF <- factor(ab$classF, levels=l, ordered=TRUE)

# list all parameter names, and split individual repseqs apart, since they will
# be in their own subdirectory; remove class and eigenvector columns
l <- colnames(ab)
l2 <- l[grep("repseq\\.only", l)]
l <- l[grep("^(NA_|[AB]_|class|HUVEC|IMR90|huvec\\.repseq\\.only)", l, invert=TRUE)]
# generate filenames for all parameters
f <- c(paste0("plot/", l, ".pdf"), paste0("plot/repseq/", l2, ".pdf"))
l <- c(l, l2)
# only generate missing plots
i <- which(file.access(f, 4) != 0)
l <- l[i]
f <- f[i]
# generate two data.frames based on which eigenvector we use, then split into
# sub-data.frames for each parameter; one list will be for raw eigenvector, and
# one for eigenvector without flanking regions (as defined in tabs.R); we will
# want to generate plots for both for each parameter; foreach will then take
# both lists and manage distributing them to workers
# for both lists, remove values in NA_ classes (erroneous data), and generate
# anyA/anyB and alwaysA/alwaysB columns for scatterplots
abf <- ab %>%
	filter(!grepl("^NA_", class)) %>%
	mutate(class=droplevels(class),
		class2=ifelse(grepl("^A_", class), "A", "B"),
		class3=ifelse(grepl("^A_.*AlwaysA", class), "A",
			ifelse(grepl("^B_.*AlwaysB", class), "B",
			NA)),
		classc=scattersplit(class))
abf <- lapply(l, function(x){
	abf %>%
		select(class, class2, class3, classc, HUVEC, !!sym(x))
})
abnf <- ab %>%
	filter(!grepl("^NA_", classF)) %>%
	mutate(classF=droplevels(classF),
		class2=ifelse(grepl("^A_", classF), "A", "B"),
		class3=ifelse(grepl("^A_.*AlwaysA", classF), "A",
			ifelse(grepl("^B_.*AlwaysB", classF), "B",
			NA)),
		classc=scattersplit(classF))
abnf <- lapply(l, function(x){
	abnf %>%
		select(classF, class2, class3, classc, HUVECnoflank, !!sym(x))
})
nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
l <- foreach(i=l, f=f, abf=abf, abnf=abnf, .inorder=FALSE, .multicombine=TRUE, .packages=c("ggplot2", "ggridges", "grid", "gridExtra", "dplyr")) %dopar% {
	theme_set(th)
	g1 <- ggviolin(abf, i, "class")
	g2 <- ggviolin(abnf, i, "classF")
	g3 <- ggridge(abf, i, "class")
	g4 <- ggridge(abnf, i, "classF")
	g5 <- ggscatter(abf, i, "HUVEC")
	g6 <- ggscatter(abnf, i, "HUVECnoflank")
	g <- grid.arrange(g1, g2, g3, g4, g5, g6)
	ggsave(f, g, width=24, height=20)
	g7 <- ggdensity2d(abf, i, "HUVEC")
	if(!is.null(g7)){
		u <- lapply(levels(abf$classc), function(x) ggdensity2dxor(abf, i, "HUVEC", x))
		g <- grid.arrange(g7, u[[1]], u[[2]], u[[3]], u[[4]], u[[5]], u[[6]], u[[7]],
			u[[8]], u[[9]], u[[10]], u[[11]])
		ggsave(sub("pdf$", "exploded.pdf", f), g, width=24, height=20)
	}
	g8 <- ggdensity2d(abnf, i, "HUVECnoflank")
	if(!is.null(g8)){
		v <- lapply(levels(abnf$classc), function(x) ggdensity2dxor(abnf, i, "HUVECnoflank", x))
		g <- grid.arrange(g8, v[[1]], v[[2]], v[[3]], v[[4]], v[[5]], v[[6]], v[[7]],
			v[[8]], v[[9]], v[[10]], v[[11]])
		ggsave(sub("pdf$", "nfexploded.pdf", f), g, width=24, height=20)
	}
}

# generate model plots
ab <- read.table("tabs/lmparms.tsv.gz", header=TRUE)
l <- list.dirs("lm", full.names=FALSE)[-1]
f <- paste0("lm/", l, "/plot.pdf")
# only generate missing plots
i <- which(file.access(f, 4) != 0)
l <- l[i]
f <- f[i]
# slice ab into a list of data.frame with weighted parameters and model fit and error
coef <- lapply(l, function(x) read.table(paste0("lm/", x, "/coef.tsv"), header=TRUE))
abl <- lapply(coef, function(x){
	x <- x[,1:2]
	d <- select(ab, eigenvectornf, !!!syms(as.character(x[-1,1])))
	for(i in 2:ncol(x))
		d[,i] <- d[,i] * x[i,2]
	d$fitted <- x[1,2] + rowSums(select(d, -1))
	d$res <- d$eigenvectornf - d$fitted
	d
})
# generate plots
l <- foreach(ab=abl, f=f, .inorder=FALSE, .multicombine=TRUE, .packages=c("ggplot2", "grid", "gridExtra", "dplyr")) %dopar% {
	theme_set(th)
	g1 <- ggfit(ab, "res", "Model residuals versus fitted values")
	g2 <- ggfit(ab, "eigenvectornf",
		"Observed eigenvector values versus fitted values", method="lm")
	g3 <- ggqnorm(ab)
	write.pdf(grid.arrange(g1, g2, g3), f, width=24, height=20)
}
stopCluster(cl)
