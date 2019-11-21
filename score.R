library(dplyr)

press <- function(lm){
	pred.res <- residuals(lm) / (1 - lm.influence(lm)$hat)
	sum(pred.res ^ 2)
}

predrsq <- function(lm){
	1 - press(lm) / sum(anova(lm)$"Sum Sq")
}

model <- function(ab, dir, expl){
	dir.create(dir)
	fx <- quote(paste0("HUVECnoflank ~", paste0(expl, collapse="+")))
	m <- lm(eval(fx), ab)
	sink(paste0(dir, "summary.txt"))
	print(summary(m))
	sink()
	pdf(paste0(dir, "plots.pdf"), width=12.1, height=9.7)
	layout(matrix(c(1,2,3,4),2,2))
	plot(m)
	dev.off()
	gfd <- gzfile(paste0(dir, "noflank.bedgraph.gz"), "wb")
	cat(paste0("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200 name=", dir, "noflank\n"),
		file=gfd)
	ab %>%
		filter(!is.na(HUVECnoflank)) %>%
		select(chr, start, end) %>%
		mutate(v=m$fitted.values) %>%
		write.table(gfd,
			col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	close(gfd)
	gfd <- gzfile(paste0(dir, "full.bedgraph.gz"), "wb")
	cat(paste0("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200 name=", dir, "full\n"),
		file=gfd)
	e <- sapply(seq_along(expl), function(i) ab[,expl[i]] * m$coefficients[1+i])
	ab %>%
		mutate(v=apply(e, 1, function(x) m$coefficients[1] + sum(x))) %>%
		filter(!is.na(HUVEC)) %>%
		select(chr, start, end, v) %>%
		write.table(gfd,
			col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	close(gfd)
	m
}

ab <- read.table("tabs/counts.tsv.gz", header=TRUE)

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
		"huvec.groseq.score"
	), list(
		"huvec.h3k27ac",
		"huvec.h3k4me3",
		"huvec.dhs.rep1",
		"huvec.silencer",
		"huvec.groseq.score",
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
	)
)

l <- apply(expand.grid(epiparms, seqparms), 1, unlist)
f <- apply(expand.grid(seq_along(epiparms)-1, seq_along(seqparms)-1), 1, function(x){
	paste0("score/m", x[1], ".", x[2], "/")
})
l <- lapply(seq_along(l), function(i){
	model(ab, f[i], l[i])
})
if(file.access("score/summary.txt") != 0){
	lapply(l, function(x){
		data.frame(rsq=round(summary(x)$r.squared, 4),
			rsqadj=round(summary(x)$adj.r.squared, 4),
			predrsq=round(predrsq(x), 4),
			vars=paste(names(x$coefficients)[-1], collapse=", "))
	}) %>%
		bind_rows %>%
		mutate(model=paste0("m", row_number())) %>%
		select(model, everything()) %>%
		arrange(desc(rsqadj)) %>%
		write.table("score/summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
}
