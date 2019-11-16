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

# generate a/b classes
ab <- read.table("prep/ab.bed", header=TRUE) %>%
	mutate(ngene=addcol(chr, "cnt/hg19.refseq.bed.gz"),
		nact=addcol(chr, "cnt/huvec.active.promoters.bed.gz"),
		class1=ifelse(HUVEC < 0, "B", "A"),
		class2=ifelse(ngene >= 4, "highgenedensity", ifelse(ngene > 0, "normalgenedensity", "nogene")),
		class2=factor(class2, levels=c("highgenedensity", "normalgenedensity", "nogene")),
		class3=ifelse(nact > 0, "hasactive", "noactive"),
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
write.gzip(ab, "tabs/class.tsv.gz", TRUE)

for(i in list.files("cnt", pattern="*.gz", full.names=TRUE)){
	if(i %in% c("hg19.refseq.bed.gz", "huvec.active.promoters.bed.gz"))
		next
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, i))
}
write.gzip(ab, "tabs/counts.tsv.gz", TRUE)
for(i in list.files("cnt/repseq", pattern="*.gz", full.names=TRUE)){
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/repseq/", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, i))
}
write.gzip(ab, "tabs/aball.tsv.gz", TRUE)

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
