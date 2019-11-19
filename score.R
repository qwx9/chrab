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

chromcounts <- list(
	"huvec.chromhmm.1.active_promoter",
	"huvec.chromhmm.2.weak_promoter",
	"huvec.chromhmm.4.strong_enhancer",
	"huvec.chromhmm.5.strong_enhancer",
	"huvec.chromhmm.6.weak_enhancer",
	"huvec.chromhmm.7.weak_enhancer",
	"huvec.chromhmm.3.poised_promoter",
	"huvec.silencer"
)
chromfeeling <- list(
	"huvec.h3k27ac",
	"huvec.h3k4me3",
	"huvec.dhs.rep1"
)
additional <- list(
	"huvec.groseq.score",
	"huvec.repseq.proa"
)
brepseq <- list(
	"huvec.repseq.prob",
	"huvec.repseq.prob.ltrline",
	"huvec.repseq.prob.ltr",
	list(
		"huvec.repseq.psil.prob",
		"huvec.repseq.prob.ltrline"
	)
)

l <- list(
	chromcounts,
	c(chromcounts, brepseq[[1]]),
	c(chromcounts, brepseq[[2]]),
	c(chromcounts, brepseq[[3]]),
	c(chromcounts, brepseq[[4]]),
	c(chromcounts, brepseq[[1]], additional[[1]]),
	c(chromcounts, brepseq[[2]], additional[[1]]),
	c(chromcounts, brepseq[[3]], additional[[1]]),
	c(chromcounts, brepseq[[4]], additional[[1]]),
	c(chromcounts, brepseq[[1]], additional[[2]]),
	c(chromcounts, brepseq[[2]], additional[[2]]),
	c(chromcounts, brepseq[[3]], additional[[2]]),
	c(chromcounts, brepseq[[4]], additional[[2]]),
	c(chromcounts, brepseq[[1]], additional),
	c(chromcounts, brepseq[[2]], additional),
	c(chromcounts, brepseq[[3]], additional),
	c(chromcounts, brepseq[[4]], additional),
	chromfeeling,
	c(chromfeeling, brepseq[[1]]),
	c(chromfeeling, brepseq[[2]]),
	c(chromfeeling, brepseq[[3]]),
	c(chromfeeling, brepseq[[4]]),
	c(chromfeeling, brepseq[[1]], additional[[1]]),
	c(chromfeeling, brepseq[[2]], additional[[1]]),
	c(chromfeeling, brepseq[[3]], additional[[1]]),
	c(chromfeeling, brepseq[[4]], additional[[1]]),
	c(chromfeeling, brepseq[[1]], additional[[2]]),
	c(chromfeeling, brepseq[[2]], additional[[2]]),
	c(chromfeeling, brepseq[[3]], additional[[2]]),
	c(chromfeeling, brepseq[[4]], additional[[2]]),
	c(chromfeeling, brepseq[[1]], additional),
	c(chromfeeling, brepseq[[2]], additional),
	c(chromfeeling, brepseq[[3]], additional),
	c(chromfeeling, brepseq[[4]], additional)
)
l <- lapply(l, unlist)
l <- l[which(sapply(seq_along(l), function(i) file.access(paste0("score/m", i)) != 0))]
l <- lapply(seq_along(l), function(i){
	model(ab, paste0("score/m", i, "/"), unlist(l[[i]]))
})
if(file.access("score/summary.txt") != 0){
	lapply(l, function(x){
		data.frame(rsq=round(summary(x)$r.squared, 2),
			rsqadj=round(summary(x)$adj.r.squared, 2),
			predrsq=round(predrsq(x), 2),
			vars=paste(names(x$coefficients)[-1], collapse=", "))
	}) %>%
		bind_rows %>%
		mutate(model=paste0("m", row_number())) %>%
		select(model, everything()) %>%
		arrange(desc(rsqadj)) %>%
		write.table("score/summary.txt", sep="\t", quote=FALSE, row.names=FALSE)
}
