# generate linear models from counted parameters, diagnostic plots and summaries for each
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)

# predicted r squared calculation
press <- function(lm){
	pred.res <- residuals(lm) / (1 - lm.influence(lm)$hat)
	sum(pred.res ^ 2)
}
predrsq <- function(lm){
	1 - press(lm) / sum(anova(lm)$"Sum Sq")
}

# generate model and model plots; takes the model data.frame, output directory
# and model name
# model is created on eigenvector with flanking regions filtered out, and is
# then run against raw eigenvector; bedgraphs for both are made
model <- function(ab, dir, name){
	# create model directory
	dir.create(dir)
	# list parameter names
	expl <- colnames(ab)[-c(1:5)]
	# generate linear model formula
	fx <- quote(paste0("HUVECnoflank ~", paste0(expl, collapse="+")))
	# compute model
	m <- lm(eval(fx), ab)

	# write out lm summary
	sink(paste0(dir, "summary.txt"))
	print(summary(m))
	sink()

	# get model diagnostic plots
	pdf(paste0(dir, "plots.pdf"), width=12.1, height=9.7)
	layout(matrix(c(1,2,3,4),2,2))
	plot(m)
	dev.off()

	# attempt for a better residuals vs fitted plot
	pdf(paste0(dir, "res.vs.fitted.pdf"), width=12.1, height=9.7)
	g <- data.frame(Fitted=m$fitted.values, Residuals=m$residuals) %>%
		ggplot(aes(Fitted, Residuals)) +
			geom_point(na.rm=TRUE, color="darkblue", alpha=0.2, shape=20) +
			geom_smooth(se=FALSE, color="red", size=0.4) +
			ggtitle("Model residuals versus fitted values")
	print(g)
	dev.off()

	# attempt for a better qqplot
	pdf(paste0(dir, "qqnorm.pdf"), width=12.1, height=9.7)
	g <- data.frame(Residuals=m$residuals) %>%
		ggplot(aes(sample=Residuals)) +
			stat_qq_line() +
			stat_qq(color="darkblue", size=0.7) +
			ggtitle("Quantile-quantile plot of residuals") +
			xlab("Theoretical normal distribution quantiles") +
			ylab("Standardized residuals")
	print(g)
	dev.off()

	# predict A/B profile on raw eigenvector; must use same variable names
	# if using predict()
	pred <- ab %>%
		mutate(HUVECnoflank=HUVEC) %>%
		select(-HUVEC)
	v <- predict(m, pred)

	# attempt at a observed vs predicted plot
	pdf(paste0(dir, "obs.vs.fitted.pdf"), width=12.1, height=9.7)
	g <- pred %>%
		rename(Observed=HUVECnoflank) %>%
		mutate(Fitted=v) %>%
		ggplot(aes(Fitted, Observed)) +
			geom_point(na.rm=TRUE, color="darkblue", alpha=0.2, shape=20) +
			geom_smooth(method="lm", na.rm=TRUE, se=FALSE, color="red", size=0.4) +
			ggtitle(paste0("Observed eigenvector values versus values predicted by model ",
				name, " (R²adj=", round(summary(m)$adj.r.squared, 2), ")"))
	print(g)
	dev.off()

	# make bedgraph with prediction on eigenvector without flanking regions
	gfd <- gzfile(paste0(dir, "noflank.bedgraph.gz"), "wb")
	# prepend header
	cat(paste0("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200 name=", dir, "noflank\n"),
		file=gfd)
	# add predicted eigenvector value for bins without missing data
	ab %>%
		filter(!is.na(HUVECnoflank)) %>%
		select(chr, start, end) %>%
		mutate(v=m$fitted.values) %>%
		write.table(gfd,
			col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	close(gfd)

	# make bedgraph with prediction on raw eigenvector
	gfd <- gzfile(paste0(dir, "full.bedgraph.gz"), "wb")
	cat(paste0("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200 name=", dir, "full\n"),
		file=gfd)
	ab %>%
		mutate(v=v) %>%
		filter(!is.na(HUVEC)) %>%
		select(chr, start, end, v) %>%
		write.table(gfd,
			col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	close(gfd)

	# return model object for further processing
	m
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
	pdf("score/cor.pdf", width=12.1, height=9.7)
	# convert matrix to data.frame, then use gather to make a tidy
	# data.frame suitable for ggplot2, then draw correlation heatmap
	g <- c %>%
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
	print(g)
	dev.off()
}

# epigenomic and genomic parameter lists to cross-combine
epiparms <- list(
	list(
		NULL
	), list(
		"huvec.chromhmm.merged.strong.promoters",
		"huvec.chromhmm.merged.weak.promoters",
		"huvec.chromhmm.any.strong.enhancers"
	), list(
		"huvec.chromhmm.merged.strong.promoters",
		"huvec.chromhmm.merged.weak.promoters",
		"huvec.chromhmm.any.strong.enhancers",
		"huvec.chromhmm.3.poised_promoter",
		"huvec.silencer"
	), list(
		"huvec.chromhmm.merged.strong.promoters",
		"huvec.chromhmm.merged.weak.promoters",
		"huvec.chromhmm.any.strong.enhancers",
		"huvec.dhs.rep1"
	), list(
		"huvec.chromhmm.merged.strong.promoters",
		"huvec.chromhmm.merged.weak.promoters",
		"huvec.chromhmm.any.strong.enhancers",
		"huvec.chromhmm.3.poised_promoter",
		"huvec.silencer",
		"huvec.dhs.rep1"
	), list(
		"huvec.h3k27ac",
		"huvec.h3k4me3",
		"huvec.dhs.rep1",
		"huvec.silencer",
		"huvec.groseq.capped"
	), list(
		"huvec.h3k27ac",
		"huvec.h3k4me3",
		"huvec.dhs.rep1",
		"huvec.silencer",
		"huvec.groseq.capped",
		"huvec.chromhmm.3.poised_promoter"
	), list(
		"huvec.dhs.rep1"
	)
)
seqparms <- list(
	list(
		NULL
	), list(
		"huvec.repseq.prob"
	), list(
		"huvec.repseq.prob.te"
	), list(
		"huvec.repseq.prob.ltr",
		"huvec.repseq.prob.l1long",
		"huvec.repseq.psil.prob"
	), list(
		"huvec.repseq.prob.ltr"
	), list(
		"huvec.repseq.proa",
		"huvec.repseq.prob"
	), list(
		"huvec.repseq.proa.te",
		"huvec.repseq.prob.te"
	), list(
		"huvec.repseq.proa.te",
		"huvec.repseq.prob.ltr",
		"huvec.repseq.prob.l1long",
		"huvec.repseq.psil.prob"
	), list(
		"huvec.repseq.proa.te",
		"huvec.repseq.prob.ltr"
	), list(
		"hg19.refseq"
	), list(
		"hg19.gc"
	), list(
		"huvec.repseq.proa"
	), list(
		"huvec.repseq.proa",
		"hg19.refseq"
	), list(
		"huvec.repseq.proa",
		"huvec.repseq.prob",
		"hg19.refseq"
	), list(
		"hg19.refseq",
		"hg19.gc"
	), list(
		"huvec.repseq.proa",
		"huvec.repseq.prob",
		"hg19.gc"
	), list(
		"huvec.repseq.proa.nonte",
		"huvec.repseq.prob.sr.lc.sat"
	), list(
		"huvec.repseq.proa.nonte",
		"huvec.repseq.prob"
	), list(
		"huvec.repseq.proa",
		"huvec.repseq.prob.sr.lc.sat"
	), list(
		"huvec.repseq.proa.te",
		"huvec.repseq.proa.nonte",
		"huvec.repseq.prob.te",
		"huvec.repseq.prob.sr.lc.sat"
	)
)

# read table of all counts, only keep useful columns
ab <- read.table("tabs/counts.tsv.gz", header=TRUE) %>%
	select(chr, start, end, HUVEC, HUVECnoflank, !!!syms(unique(unlist(c(epiparms, seqparms)))))
# normalize parameter columns range to [0;1] (doesn't affect prediction
# efficiency, but coeffs are more easily interpretable)
ab[,6:ncol(ab)] <- apply(ab[,6:ncol(ab)], 2, function(x) x / max(x, na.rm=TRUE))

# export csv with eigenvector and transformed parameters for neural network
ab %>%
	select(-chr, -start, -end, -HUVECnoflank) %>%
	rename(eigenvector=HUVEC) %>%
	write.csv(file="score/params.csv", quote=FALSE, row.names=FALSE)

# generate correlation heatmap for all parameters
ggcor(ab[,-c(1:3,5)])

# get list of all possible combinations between the two lists, and names and
# filenames of these combinations
l <- apply(expand.grid(epiparms, seqparms), 1, unlist)
n <- apply(expand.grid(seq_along(epiparms)-1, seq_along(seqparms)-1), 1, function(x){
	paste0("m", x[1], ".", x[2])
})
f <- paste0("score/", n, "/")
# remove first combination (NULL, NULL)
l <- l[-1]
n <- n[-1]
f <- f[-1]
# generate a list of table slices for each combination, which foreach will
# distribute among workers
abl <- lapply(l, function(x){
	select(ab, chr, start, end, HUVEC, HUVECnoflank, !!!syms(as.character(x)))
})

cl <- makeCluster(detectCores())
registerDoParallel(cl)
l <- foreach(ab=abl, f=f, n=n, .multicombine=TRUE, .inorder=FALSE, .packages=c("ggplot2", "dplyr")) %dopar% {
	# generate model and make plots; return a list of model metrics and parameters
	m <- model(ab, f, n)
	list(data.frame(model=as.character(n),
		rsq=round(summary(m)$r.squared, 4),
		rsqadj=round(summary(m)$adj.r.squared, 4),
		predrsq=round(predrsq(m), 4),
		vars=paste(names(m$coefficients)[-1], collapse=", "),
		stringsAsFactors=FALSE),
		data.frame(model=n, as.list(m$coefficients[-1]),
		stringsAsFactors=FALSE))
}
stopCluster(cl)

# write a summary of each model's performance, sorted by adjusted rsquare
lapply(l, function(x) x[[1]]) %>%
	bind_rows %>%
	arrange(desc(rsqadj)) %>%
	write.table("score/summary.txt", sep="\t", quote=FALSE, row.names=FALSE)

# empty template data.frame for a table of parameters for each model
abz <- ab %>%
	mutate(model="") %>%
	select(model, everything()) %>%
	select(-chr, -start, -end, -HUVEC, -HUVECnoflank) %>%
	filter(rep(FALSE))
# export a table for parameter used across all models
lapply(l, function(x) full_join(abz, x[[2]], by=colnames(x[[2]]))) %>%
	bind_rows %>%
	write.table("score/params.txt", sep="\t", quote=FALSE, row.names=FALSE)
