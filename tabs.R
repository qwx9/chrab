require(dplyr)

addcol <- function(chr, f){
	read.table(f) %>%
		filter(V1 %in% unique(chr)) %>%
		pull(V4)
}

# generate a/b classes
ab <- read.table("prep/ab.bed", header=TRUE) %>%
	mutate(ngene=addcol(chr, "cnt/hg19w.hg19.refseq.bed.gz")) %>%
	mutate(nact=addcol(chr, "cnt/hg19w.huvec.proa.genes.bed.gz")) %>%
	mutate(class1=ifelse(HUVEC < 0, "B", "A"),
		class2=ifelse(ngene >= 4, "highgenedensity", ifelse(ngene > 0, "normalgenedensity", "nogene")),
		class3=ifelse(nact > 0, "hasactive", "noactive"),
		class=paste(class1, class2, AorBvec, class3, sep="_")) %>%
	select(-class1, -class2, -class3)
write.table(ab, "tabs/class.tsv", sep="\t", row.names=FALSE, quote=FALSE)

for(i in list.files("cnt")){
	if(i %in% c("hg19w.hg19.refseq.bed.gz", "hg19w.huvec.proa.genes.bed.gz"))
		next
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/[^\\.]+\\.", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, paste0("cnt/", i)))
}
gfd <- gzfile("tabs/aball.tsv.gz", "w")
write.table(ab, gfd, sep="\t", row.names=FALSE, quote=FALSE)
close(gfd)
