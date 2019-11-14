require(dplyr)
source("lib.R")

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
		class=paste(class1, class2, AorBvec, class3, sep="_"))

as.data.frame(table(ab$class1, useNA="always")) %>%
	select(type=Var1, count=Freq) %>%
	arrange(type) %>%
	filter(count != 0) %>%
	write.table("tabs/bins.by1.tsv", sep="\t", row.names=FALSE, quote=FALSE)
as.data.frame(table(ab$class1, ab$class2, useNA="always")) %>%
	select(type=Var1, genedensity=Var2, count=Freq) %>%
	mutate(genedensity=factor(genedensity, levels=levels(genedensity)[c(1,3,2)])) %>%
	arrange(type, genedensity) %>%
	filter(count != 0) %>%
	write.table("tabs/bins.by2.tsv", sep="\t", row.names=FALSE, quote=FALSE)
as.data.frame(table(ab$class1, ab$class2, ab$AorBvec, useNA="always")) %>%
	select(type=Var1, genedensity=Var2, pc1=Var3, count=Freq) %>%
	mutate(genedensity=factor(genedensity, levels=levels(genedensity)[c(1,3,2)])) %>%
	arrange(type, genedensity, pc1) %>%
	filter(count != 0) %>%
	write.table("tabs/bins.by3.tsv", sep="\t", row.names=FALSE, quote=FALSE)
as.data.frame(table(ab$class1, ab$class2, ab$AorBvec, ab$class3, useNA="always")) %>%
	select(type=Var1, genedensity=Var2, pc1=Var3, active=Var4, count=Freq) %>%
	mutate(genedensity=factor(genedensity, levels=levels(genedensity)[c(1,3,2)])) %>%
	arrange(type, genedensity, pc1, active) %>%
	filter(count != 0) %>%
	write.table("tabs/bins.by4.tsv", sep="\t", row.names=FALSE, quote=FALSE)

ab <- ab %>%
	select(-class1, -class2, -class3)
write.table(ab, "tabs/class.tsv", sep="\t", row.names=FALSE, quote=FALSE)

for(i in list.files("cnt")){
	if(i %in% c("hg19w.hg19.refseq.bed.gz", "hg19w.huvec.proa.genes.bed.gz"))
		next
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/[^\\.]+\\.", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, paste0("cnt/", i)))
}
write.gzip(ab, "tabs/aball.tsv.gz")
ab %>%
	select(-contains("repseq.only")) %>%
	write.table("tabs/counts.tsv", sep="\t", row.names=FALSE, quote=FALSE)

for(i in unique(ab$class)){
	ab %>%
		mutate(v=ifelse(class==i, 1, 0)) %>%
		select(chr, start, end, v) %>%
		write.gzip(paste0("cnt/hg19w.", i, ".bed.gz"))
}
