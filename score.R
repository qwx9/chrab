library(dplyr)

model <- function(ab, pred, expl){
	fx <- quote(paste0(pred, "~", paste0(expl, collapse="+")))
	m <- lm(eval(fx), ab)
	f <- paste0("score/", pred, ".", paste0(expl, collapse="."))
	f <- gsub("(huvec|chromhmm\\.[0-9]|repseq)\\.", "", f)
	if(file.access(paste0(f, ".txt")) == 0)
		return(NULL)
	sink(paste0(f, ".txt"))
	print(summary(m))
	sink()
	pdf(paste0(f, ".pdf"), width=12.1, height=9.7)
	layout(matrix(c(1,2,3,4),2,2))
	plot(m)
	dev.off()
	gfd <- gzfile(paste0(f, ".bedgraph.gz"), "wb")
	cat(paste0("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200 name=", f, "\n"),
		file=gfd)
	ab %>%
		filter(!is.na(!!sym(pred))) %>%
		select(chr, start, end) %>%
		mutate(v=m$fitted.values) %>%
		write.table(gfd,
			col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	close(gfd)
}

ab <- read.table("tabs/counts.tsv.gz", header=TRUE)

proa <- list(
	list(
		"huvec.chromhmm.1.active_promoter",
		"huvec.chromhmm.2.weak_promoter",
		"huvec.chromhmm.3.poised_promoter",
		"huvec.chromhmm.4.strong_enhancer",
		"huvec.chromhmm.5.strong_enhancer",
		"huvec.chromhmm.6.weak_enhancer",
		"huvec.chromhmm.7.weak_enhancer",
		"huvec.silencer"
	),
	list(
		"huvec.h3k27ac",
		"huvec.h3k4me3",
		"huvec.dhs.rep1"
	)
)
prox <- list(
	"hg19.gc",
	"huvec.groseq.score",
	"huvec.repseq.proa"
)
prob <- list(
	"huvec.repseq.ltrline",
	"huvec.repseq.ltr",
	"huvec.repseq.l1",
	"huvec.repseq.l1.long",
	"huvec.repseq.prob",
	"huvec.repseq.prob.ltrline",
	"huvec.repseq.prob.ltr",
	"huvec.repseq.prob.line",
	"huvec.repseq.prob.sr.lc.sat",
	"huvec.repseq.psil.anyb",
	"huvec.repseq.psil.prob",
	"huvec.repseq.psil.subb"
)
l <- lapply(proa, function(x) lapply(prob, function(y){
	model(ab, "HUVEC", c(x, y))
	model(ab, "HUVECnoflank", c(x, y))
	lapply(prox, function(z){
		model(ab, "HUVEC", c(x, y, z))
		model(ab, "HUVECnoflank", c(x, y, z))
	})
}))
