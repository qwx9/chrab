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

mkrefseq <- function(f){
	read.table("hg19/ncbiRefSeq.txt.gz") %>%
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
		mutate(score=rep(0)) %>%
		select(chr, start, end, gene, score, strand) %>%
		write.gzip(f)
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
			cnt[n] <- as.integer(i[2]) / 20
			nt[n] <- 5
		}else if(as.integer(i[1]) - 1 > end[n]){
			n <- n + 1
			chr[n] <- chr[n-1]
			start[n] <- start[n-1] + 100000
			end[n] <- end[n-1] + 100000
			cnt[n] <- as.integer(i[2]) / 20
			nt[n] <- 5
		}else{
			cnt[n] <- cnt[n] + as.integer(i[2]) / 20
			nt[n] <- nt[n] + 5
		}
	}
	list(chr=chr[1:n], start=start[1:n], end=end[1:n], n=cnt[1:n], nt=nt[1:n], sq=rep(sq, n))
}

mkgc5 <- function(f){
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
	chrom <- read.table("prep/hg19.txt", header=TRUE)
	df <- df %>%
		mutate(chr=as.character(l[chr])) %>%
		filter(chr %in% chrom$chrom) %>%
		group_by(chr, start, end) %>%
		summarize(n=sum(n), nt=sum(nt)) %>%
		ungroup
	ends <- df %>%
		group_by(chr) %>%
		summarize(len=max(end))
	l <- unlist(lapply(chrom$chr, function(x) if(chrom[chrom$chrom==x,2] > ends[ends$chr==x,2]) as.character(x)))
	df2 <- bind_rows(lapply(l, function(i){
		a <- chrom[chrom$chrom == i,]
		b <- ends[ends$chr == i,]
		b <- data.frame(chr=factor(i, levels=levels(chrom$chrom)), start=b$len, end=b$len+100000, n=NA, nt=NA)
		r <- b
		while(r$end < a$size){
			r$start <- r$end
			r$end <- r$end + 100000
			b <- rbind(b, r)
		}
		b
	}))
	df %>%
		rbind(df2) %>%
		mutate(gc=n/nt, chr0=ifelse(nchar(chr) == 4 & substr(chr, 4, 4) %in% 0:9, paste0("chr0", substr(chr, 4, 4)), chr)) %>%
		arrange(chr0, start) %>%
		select(chr, start, end, gc, n, nt, -chr0) %>%
		write.gzip(f)
}

mkchr <- function(f){
	o <- src_mysql("hg19", "genome-euro-mysql.soe.ucsc.edu", username="genome")
	t <- o %>%
		tbl("chromInfo") %>%
		select(chrom, size) %>%
		collect %>%
		filter(grepl("chr[1-9XY][0-9]?$", chrom)) %>%
		mutate(chr=ifelse(nchar(chrom) == 4 & substr(chrom, 4, 4) %in% 0:9, paste0("chr0", substr(chrom, 4, 4)), chrom)) %>%
		arrange(chr) %>%
		select(-chr) %>%
		write.table(f, sep="\t", row.names=FALSE, quote=FALSE)
}

mkab <- function(f){
	read_xlsx("gf/Table-AouBouAlways.xlsx", col_types="text") %>%
		select(chr, start, end, AorBvec, HUVEC, IMR90) %>%
		mutate(start=format(as.integer(start)-1, scientific=FALSE, trim=TRUE)) %>%
		write.table(f, sep="\t", quote=FALSE, row.names=FALSE)
}

mkhuvecrep <- function(f){
	read_xlsx("gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx", col_types="text",
		.name_repair=~gsub(" ", "_", .x)) %>%
		select(-Rank_CorrB) %>%
		write.table(f, sep="\t", quote=FALSE, row.names=FALSE)
}

mkpromenh <- function(f){
	read.table("hg19/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz", fill=TRUE) %>%
		filter(V3 %in% c("promoter", "enhancer")) %>%
		select(V1, V4, V5, V9) %>%
		mutate(V1=paste0("chr", V1), V9=strsplit(as.character(V9), "=|;")[[1]][2]) %>%
		write.gzip(f)
}

mkprom <- function(f){
	read.table("hg19/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz", fill=TRUE) %>%
		filter(V3 == "promoter") %>%
		select(V1, V4, V5, V9) %>%
		mutate(V1=paste0("chr", V1), V9=strsplit(as.character(V9), "=|;")[[1]][2]) %>%
		write.gzip(f)
}

mkgroseq <- function(f){
	lapply(c("huvec/groseq.rep1.bedGraph.gz",
		"huvec/groseq.rep2.bedGraph.gz",
		"huvec/groseq.rep3.bedGraph.gz",
		"huvec/groseq.rep4.bedGraph.gz"),
		function(x){
			x <- read.table(x, skip=1, sep="\t", fill=TRUE)
			x <- subset(x, V1 != "chrM")
			i <- grep("^track ", x$V1)
			x <- droplevels(x[-i,])
			y <- x[i:nrow(x),]
			y$strand <- "-"
			x <- x[2:(i-1),]
			x$strand <- "+"
			list(x, y)
		}) %>%
		unlist(recursive=FALSE) %>%
		bind_rows %>%
		mutate(name=rep("")) %>%
		select(V1, V2, V3, name, V4, strand) %>%
		arrange(V1, V2) %>%
		write.gzip(f)
}

mksilencer <- function(f){
	zfd <- unz("huvec/ovcharenko2019.supp.tab.s2.txt.zip", "Supplemental_Table_S2.txt")
	read.table(zfd, skip=1) %>%
		filter(row_number() %in% contains("HUVEC", vars=V1)) %>%
		mutate(pos=strsplit(as.character(V2), ":|-"),
			chr=sapply(pos, function(x) x[1]),
			start=sapply(pos, function(x) x[2]),
			end=sapply(pos, function(x) x[3])) %>%
		select(chr, start, end) %>%
		arrange(chr, as.integer(start)) %>%
		write.gzip(f)
}

mkconv <- function(){
	l <- list(
		list(f="prep/hg19.refseq.bed.gz", fn=mkrefseq),
		list(f="prep/hg19.txt", fn=mkchr),
		list(f="prep/hg19w.hg19.gc5base.bed.gz", fn=mkgc5),
		list(f="prep/ab.bed", fn=mkab),
		list(f="prep/huvec.repseq.tsv", fn=mkhuvecrep),
		list(f="prep/hg19.promenh.gff.gz", fn=mkpromenh),
		list(f="prep/hg19.prom.gff.gz", fn=mkprom),
		list(f="prep/huvec.groseq.allrep.bed.gz", fn=mkgroseq),
		list(f="prep/huvec.silencer.bed.gz", fn=mksilencer)
	)
	for(i in l)
		if(file.access(i$f, 4) != 0)
			i$fn(i$f)
}

mkchromhmm <- function(){
	l <- c("1_Active_Promoter",
		"2_Weak_Promoter",
		"3_Poised_Promoter",
		"4_Strong_Enhancer",
		"5_Strong_Enhancer",
		"6_Weak_Enhancer",
		"7_Weak_Enhancer")
	f <- paste0("prep/huvec.chromhmm.", tolower(sub("_", ".", sub("[0-9]+_", "", l))), ".bed.gz")
	i <- which(file.access(f, 4) != 0)
	if(length(i) != 0){
		x <- read.table("hg19/wgEncodeBroadHmmHuvecHMM.bed.gz")
		x <- sapply(i, function(i){
			x %>%
				filter(V4 == l[i]) %>%
				select(V1, V2, V3) %>%
				write.gzip(f[i])
			NULL
		})
	}
}

mkrepseq <- function(){
	repmask <- read.table("hg19/hg19.fa.out.gz", skip=2) %>%
		select(chr=V5, start=V6, end=V7, id=V10)
	repcor <- read.table("prep/huvec.repseq.tsv", header=TRUE)

	l <- list(
		list(f="prep/huvec.repseq.ltrline.bed.gz", q=quote(repcor %>% filter(Class %in% c("LTR", "LINE")))),
		list(f="prep/huvec.repseq.ltr.bed.gz", q=quote(repcor %>% filter(Class == "LTR"))),
		list(f="prep/huvec.repseq.l1.bed.gz", q=quote(repcor %>% filter(Class == "LINE" & Family == "L1"))),
		list(f="prep/huvec.repseq.prob.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01))),
		list(f="prep/huvec.repseq.prob.ltrline.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class %in% c("LTR", "LINE")))),
		list(f="prep/huvec.repseq.prob.ltr.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class == "LTR"))),
		list(f="prep/huvec.repseq.prob.line.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class == "LINE"))),
		list(f="prep/huvec.repseq.prob.sr.lc.sat.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class %in% c("Simple_repeat", "Low_complexity", "Satellite")))),
		list(f="prep/huvec.repseq.proa.bed.gz", q=quote(repcor %>% filter(Moyenne_corA >= 0.04)))
	)
	f <- paste0("prep/huvec.repseq.only.",
		gsub("\\(|\\)", "", gsub("/", "-", repcor$Name)), ".bed.gz")
	f <- c(sapply(l, function(x) x$f), f)
	n <- c(sapply(l, function(x) as.character(eval(x$q)$Name)), as.character(repcor$Name))
	i <- unlist(sapply(seq_along(f), function(i) if(file.access(f[i], 4) == 0) i else NULL))
	if(!is.null(i)){
		f <- f[-i]
		n <- n[-i]
	}
	nc <- detectCores()
	cl <- makeCluster(nc)
	registerDoParallel(cl)
	l <- foreach(f=f, n=n, .inorder=FALSE, .export="write.gzip") %dopar%
		write.gzip(repmask[repmask$id %in% n,], f)
	stopCluster(cl)

	if(file.access("prep/huvec.repseq.l1.long.bed.gz", 4) != 0)
		read.table("prep/huvec.repseq.l1.bed.gz") %>%
			filter(V3-V2 > 5000) %>%
			write.gzip("prep/huvec.repseq.l1.long.bed.gz")
}

mkconv()
mkchromhmm()
mkrepseq()
