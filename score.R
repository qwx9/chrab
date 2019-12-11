library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)

press <- function(lm){
	pred.res <- residuals(lm) / (1 - lm.influence(lm)$hat)
	sum(pred.res ^ 2)
}

predrsq <- function(lm){
	1 - press(lm) / sum(anova(lm)$"Sum Sq")
}

model <- function(ab, dir){
	dir.create(dir)
	expl <- colnames(ab)[-c(1:5)]
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

ggcor <- function(ab){
	c <- cor(ab, use="complete.obs")
	hc <- hclust(as.dist((1-c)/2))
	c <- c[hc$order, hc$order]
	write.table(c, "score/cor.txt", sep="\t", quote=FALSE)
	c[lower.tri(c)] <- NA
	pdf("score/cor.pdf", width=12.1, height=9.7)
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

ab <- read.table("tabs/counts.tsv.gz", header=TRUE) %>%
	select(chr, start, end, HUVEC, HUVECnoflank, !!!syms(unique(unlist(c(epiparms, seqparms)))))
ab[,6:ncol(ab)] <- apply(ab[,6:ncol(ab)], 2, function(x) x / max(x, na.rm=TRUE))

ab %>%
	select(-chr, -start, -end, -HUVECnoflank) %>%
	rename(eigenvector=HUVEC) %>%
	write.csv(file="score/params.csv", quote=FALSE, row.names=FALSE)

ggcor(ab[,-c(1:3,5)])

l <- apply(expand.grid(epiparms, seqparms), 1, unlist)
n <- apply(expand.grid(seq_along(epiparms)-1, seq_along(seqparms)-1), 1, function(x){
	paste0("m", x[1], ".", x[2])
})
f <- paste0("score/", n, "/")
l <- l[-1]
n <- n[-1]
f <- f[-1]
abl <- lapply(l, function(x){
	select(ab, chr, start, end, HUVEC, HUVECnoflank, !!!syms(as.character(x)))
})

cl <- makeCluster(detectCores())
registerDoParallel(cl)
l <- foreach(ab=abl, f=f, n=n, .multicombine=TRUE, .inorder=FALSE, .packages=c("ggplot2", "dplyr")) %dopar% {
	m <- model(ab, f)
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

lapply(l, function(x) x[[1]]) %>%
	bind_rows %>%
	arrange(desc(rsqadj)) %>%
	write.table("score/summary.txt", sep="\t", quote=FALSE, row.names=FALSE)

abz <- ab %>%
	mutate(model="") %>%
	select(model, everything()) %>%
	select(-chr, -start, -end, -HUVEC, -HUVECnoflank) %>%
	filter(rep(FALSE))
lapply(l, function(x) full_join(abz, x[[2]], by=colnames(x[[2]]))) %>%
	bind_rows %>%
	write.table("score/params.txt", sep="\t", quote=FALSE, row.names=FALSE)
