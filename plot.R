# generate plots for all counts
library(dplyr)
library(tidyr)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggridges)
library(doParallel)
source("lib.R")

write.pdf <- function(x, file, width=12.1, height=9.7){
	pdf(file, width, height)
	print(x)
	dev.off()
}

# generate correlation heatmap of params in provided data.frame
ggcor <- function(ab){
        # calculate correlation between all parameters, removing NAs
        c <- cor(ab, use="complete.obs")
        # perform hierarchical clustering using the correlation as a
        # dissimilarity measure and order correlation matrix by clusters
        hc <- hclust(as.dist((1-c)/2))
        c <- c[hc$order, hc$order]
        # export raw matrix
        write.table(c, "score/cor.txt", sep="\t", quote=FALSE)
        # remove lower matrix triangle since it's redundant
        c[lower.tri(c)] <- NA
        # convert matrix to data.frame, then use gather to make a tidy
	# data.frame suitable for ggplot2, and return correlation heatmap
	# ggplot object
        c %>%
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
	# for non-integral "counts" like %gc, don't label violin plots
	# otherwise, annotate each violin with its number of observations
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

# ridgeline plots: takes a data.frame with class and parameter; parameter name;
# class variable to use
ggridge <- function(ab, var, class){
	ggplot(ab, aes(!!sym(var), !!sym(class), fill=!!sym(class))) +
		scale_fill_discrete(guide=FALSE) +
		geom_density_ridges2(na.rm=TRUE)
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
			scale_fill_gradientn(colors=c("black", "blue", "red", "yellow"),
				values=c(0, 1/50, 1/2, 1),
				breaks=unique(round(seq(0, cnt, length.out=12), 0))) +
			guides(fill=guide_colorbar(barheight=30))
	}
	# put enrichment summary on top of scatter plot
	grid.arrange(g, top=textGrob(s))
}

# read table of all counts, removing useless columns
ab <- read.table("tabs/aball.tsv.gz", header=TRUE) %>%
	select(-chr, -start, -end)

# linear models plots
# get epigenomic and genomic parameter lists
l <- lapply(c("score.eparm.tsv", "score.gparm.tsv"), read.parms)
# generate correlation heatmap for all parameters
if(file.access("score/cor.pdf", 4) != 0){
	ab %>%
		select(HUVEC, !!!syms(unique(unlist(l)))) %>%
		ggcor %>%
		write.pdf("score/cor.pdf")
}

# list all parameter names, and split individual repseqs apart, since they will
# be in their own subdirectory; remove class and eigenvector columns
l <- colnames(ab)
l2 <- l[grep("repseq\\.only", l)]
l <- l[grep("^(NA_|[AB]_|class|HUVEC|IMR90|huvec\\.repseq\\.only)", l, invert=TRUE)]
# generate filenames for all parameters
f <- c(paste0("plot/", l, ".pdf"), paste0("plot/repseq/", l2, ".pdf"))
l <- c(l, l2)
# get list of missing files, exit if there are none
i <- which(file.access(f, 4) != 0)
if(length(i) == 0)
	quit()
# only generate missing plots
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
	write.pdf(grid.arrange(g1, g2, g3, g4, g5, g6), f, width=24, height=20)
}
stopCluster(cl)
