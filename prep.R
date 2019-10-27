require(dplyr)

write.gzip <- function(x, file, ...){
	gfd <- gzfile(file, "wb")
	write.table(x, gfd, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", ...)
	close(gfd)
}

cntpref <- function(x, g){
	sapply(x, function(x){
		x <- as.integer(unlist(strsplit(x, ",")))
		if(length(x) < 2)
			return("")
		g <- as.character(g[x])
		x <- x[sapply(g, function(x){
			l <- ifelse(g != x, nchar(g), 0)
			for(i in 1:max(l)){
				y <- g[l > 0]
				if(all(ifelse(length(y)==0, TRUE, substr(x, 1, i) != substr(y, 1, i))))
					return(i-1)
				l <- l - 1
			}
			i
		}) > 0]
		paste0(x, collapse=",")
	})
}

read.table("hg19/ncbiRefSeqCurated.txt.gz") %>%
	select(chr=V3, start=V5, end=V6, strand=V4, gene=V13) %>%
	group_by(chr, gene, strand) %>%
	summarize(start=min(start), end=max(end)) %>%
	arrange(chr, start, end) %>%
	group_by(chr) %>%
	mutate(over=sapply(1:length(start), function(x){
		paste0(which(
			start[x] >= start & start[x] <= end |
			end[x] >= start & end[x] <= end |
			end >= start[x] & end <= end[x] |
			start >= start[x] & start <= end[x]), collapse=",")),
		overg=cntpref(over, gene)) %>%
	group_by(chr, overg, strand) %>%
	mutate(gene=ifelse(overg!="", paste0(as.character(gene[1]), "_", "comb"), as.character(gene))) %>%
	group_by(chr, gene, strand) %>%
	summarize(start=min(start), end=max(end)) %>%
	ungroup %>%
	arrange(chr, start, end, strand) %>%
	select(chr, start, end, strand, gene) %>%
	write.gzip("cnt/hg19.refseq.bed.gz")
