require(dplyr)
require(doParallel)
require(readxl)

gzreadlines <- function(file, n, fn){
	gfd <- gzfile(file, "r")
	sq <- 0
	nc <- detectCores()
	si <- splitIndices(nc * n, nc)
	l <- list()
	cl <- makeCluster(nc)
	registerDoParallel(cl)
	while(length(x <- readLines(gfd, nc * n)) != 0){
		if(length(x) < nc * n)
			si <- splitIndices(length(x), nc)
		x <- lapply(si, function(i) x[i])
		ll <- foreach(x=x, sid=seq(sq, sq+nc-1), .inorder=FALSE) %dopar%
			fn(x, sid)
		l <- c(l, ll)
		sq <- sq + nc
	}
	stopCluster(cl)
	close(gfd)
	bind_rows(l[sort.int(sapply(l, function(x) x$sq[1]), index.return=TRUE)$ix])
}

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

mkrefseq <- function(){
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
				start >= start[x] & start <= end[x]),
				collapse=",")}),
			overg=cntpref(over, gene)) %>%
		group_by(chr, overg, strand) %>%
		mutate(gene=ifelse(overg!="", paste0(as.character(gene[1]), "_", "comb"), as.character(gene))) %>%
		group_by(chr, gene, strand) %>%
		summarize(start=min(start), end=max(end)) %>%
		ungroup %>%
		arrange(chr, start, end, strand) %>%
		select(chr, start, end, strand, gene) %>%
		write.gzip("cnt/hg19.refseq.bed.gz")
}

gc5cnt <- function(x, sq){
	m <- length(x) / (100000 / 10)
	chr <- character(m)
	start <- integer(m)
	end <- integer(m)
	cnt <- integer(m)
	nt <- integer(m)
	n <- 0
	x <- strsplit(x, "\\s|=")
	for(i in x){
		if(i[1] == "variableStep"){
			n <- n + 1
			chr[n] <- i[3]
			start[n] <- 0
			end[n] <- 100000
			cnt[n] <- 0
			nt[n] <- 0
		}else if(n == 0){
			n <- n + 1
			k <- as.integer((as.integer(i[1]) - 1) / 100000) * 100000
			chr[n] <- paste0("chrsq", sq)
			start[n] <- k
			end[n] <- k + 100000
			cnt[n] <- 0
			nt[n] <- 0
		}else if(as.integer(i[1]) - 1 >= end[n]){
			n <- n + 1
			chr[n] <- chr[n-1]
			start[n] <- start[n-1] + 100000
			end[n] <- end[n-1] + 100000
			cnt[n] <- 0
			nt[n] <- 0
		}else{
			cnt[n] <- cnt[n] + as.integer(i[2]) / 20
			nt[n] <- nt[n] + 5
		}
	}
	list(chr=chr[1:n], start=start[1:n], end=end[1:n], n=cnt[1:n], nt=nt[1:n], sq=rep(sq, n))
}

mkgc5 <- function(){
	df <- gzreadlines("hg19/hg19.gc5Base.wigVarStep.gz", 10*100000, gc5cnt)
	u <- unique(df$chr)
	l <- vector(mode="list", length=length(u))
	names(l) <- u
	for(i in unique(df$chr)){
		if(grepl("chrsq", i))
			l[[i]] <- chr
		else{
			chr <- i
			l[[i]] <- i
		}
	}
	df %>%
		mutate(chr=as.character(l[chr])) %>%
		filter(chr %in% u[grep("chr[1-9XY][0-9]?$", u)]) %>%
		group_by(chr, start, end) %>%
		summarize(n=sum(n), nt=sum(nt)) %>%
		ungroup %>%
		mutate(gc=n/nt) %>%
		select(chr, start, end, gc, n, nt) %>%
		arrange(chr, start) %>%
		write.gzip("cnt/hg19.gc5base.bed.gz")
}

mkchr <- function(){
	o <- src_mysql("hg19", "genome-mysql.cse.ucsc.edu", username="genome")
	t <- o %>%
		tbl("chromInfo") %>%
		select(chrom, size) %>%
		collect %>%
		filter(grepl("chr[1-9XY][0-9]?$", chrom)) %>%
		mutate(chr=ifelse(nchar(chrom) == 4 & substr(chrom, 4, 4) %in% 0:9, paste0("chr0", substr(chrom, 4, 4)), chrom)) %>%
		arrange(chr) %>%
		select(-chr) %>%
		write.table("cnt/hg19.txt", sep="\t", row.names=FALSE, quote=FALSE)
}

mkab <- function(){
	read_xlsx("gf/Table-AouBouAlways.xlsx", col_types="text") %>%
		select(chr, start, end, AorBvec, HUVEC, IMR90) %>%
		mutate(start=format(as.integer(start)-1, scientific=FALSE, trim=TRUE)) %>%
		write.table("gf/ab.bed", sep="\t", quote=FALSE, row.names=FALSE)
}

mkhuvecrep <- function(){
	read_xlsx("gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx", col_types="text") %>%
		select(-"Rank CorrB") %>%
		write.table("gf/huvec.repseq.tsv", sep="\t", quote=FALSE, row.names=FALSE)
}

l <- list(
	list(f="cnt/hg19.refseq.bed.gz", fn=mkrefseq),
	list(f="cnt/hg19w.hg19.gc5base.bed.gz", fn=mkgc5),
	list(f="cnt/hg19.txt", fn=mkchr),
	list(f="prep/ab.bed", fn=mkab),
	list(f="prep/huvec.repseq.tsv", fn=mkhuvecrep)
)
for(i in l)
	if(file.access(i$f, 4) != 0)
		i$fn()
