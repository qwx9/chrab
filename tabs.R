require(dplyr)
source("lib.R")

write.counts <- function(x, file){
	l <- colnames(x)
	x <- sapply(seq_along(l), function(i){
		if(i > 1){
			i <- 1:i
			l <- l[i]
		}else{
			l <- "."
		}
		x %>%
			select(i) %>%
			table(useNA="always") %>%
			as.data.frame %>%
			rename(count=Freq) %>%
			arrange(!!!syms(l)) %>%
			write.table(paste0("tabs/", file, ".by", length(i), ".tsv"),
				sep="\t", row.names=FALSE, quote=FALSE)
	})
}

addcol <- function(chr, f){
	read.table(f) %>%
		filter(V1 %in% unique(chr)) %>%
		pull(V4)
}

ab <- read.table("prep/ab.bed", header=TRUE)
for(i in list.files("cnt", pattern="*.gz", full.names=TRUE)){
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, i))
}
cap <- quantile(ab$huvec.groseq.score, 0.75) + 5 * IQR(ab$huvec.groseq.score)
ab <- ab %>%
        mutate(huvec.groseq.capped=ifelse(huvec.groseq.score > cap, cap, huvec.groseq.score))
write.gzip(ab, "tabs/counts.tsv.gz", TRUE)
ab <- ab %>%
	mutate(class1=ifelse(HUVEC < 0, "B", "A"),
		class2=ifelse(hg19.refseq >= 4, "highgenedensity", ifelse(hg19.refseq > 0, "normalgenedensity", "nogene")),
		class2=factor(class2, levels=c("highgenedensity", "normalgenedensity", "nogene")),
		class3=ifelse(huvec.chromhmm.any.active.promoters > 0, "hasactive", "noactive"),
		class4=ifelse(HUVECnoflank < 0, "B", "A"),
		class=paste(class1, class2, AorBvec, class3, sep="_"),
		classF=paste(class4, class2, AorBvec, class3, sep="_"))
ab %>%
	select(type=class1, genedensity=class2, pc1=AorBvec, active=class3) %>%
	write.counts("bins")
ab %>%
	select(type=class4, genedensity=class2, pc1=AorBvec, active=class3) %>%
	write.counts("bins.noflank")
ab <- ab %>%
	select(-class1, -class2, -class3, -class4)
ab %>%
	select(chr, start, end, HUVEC, HUVECnoflank, class, classF) %>%
	write.gzip("tabs/class.tsv.gz", TRUE)

for(i in list.files("cnt/repseq", pattern="*.gz", full.names=TRUE)){
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/repseq/", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, i))
}
ab %>%
	select(-AorBvec, -IMR90) %>%
	select(chr, start, end, HUVEC, HUVECnoflank, class, classF, everything()) %>%
	write.gzip("tabs/aball.tsv.gz", TRUE)

n <- colnames(ab)
n <- n[grep("^(NA_|[AB]_|AorB|class|chr|start|end)", n, invert=TRUE)]
bind_rows(lapply(n, function(n){
	x <- ab[,n]
	y <- na.omit(x)
	data.frame(min=min(y), max=max(y), mean=mean(y), sd=sd(y), median=median(y),
		q1=quantile(y, 0.25), q3=quantile(y, 0.75), na=sum(is.na(x)))
})) %>%
	mutate(parm=n) %>%
	select(parm, everything()) %>%
	write.gzip("tabs/parms.tsv.gz", TRUE)

for(i in unique(ab$class)){
	ab %>%
		mutate(v=ifelse(class==i, 1, 0)) %>%
		select(chr, start, end, v) %>%
		write.gzip(paste0("cnt/", i, ".bed.gz"))
	ab %>%
		mutate(v=ifelse(classF==i, 1, 0)) %>%
		select(chr, start, end, v) %>%
		write.gzip(paste0("cnt/", i, ".noflank.bed.gz"))
}
