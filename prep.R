# data extraction, conversions and other transformations prior to counting
library(dplyr)
library(doParallel)
library(readxl)
source("lib.R")

# read compressed wig count file in parallel
# file: input file, n: chunk size, fn: chunk transform function
gzreadlines <- function(file, n, fn){
	gfd <- gzfile(file, "r")
	# prep for parallel processing: each returned chunk will have a
	# sequence number to reorder chunks and not lose any information; n
	# being a reasonable value, each worker will have a chunk of size of n
	# to process (n should be as large as possible, capped by total memory
	# usage)
	sq <- 0
	nc <- detectCores()
	# indices to split read for each worker in a balanced way
	si <- splitIndices(nc * n, nc)
	# storage for list of processed chunks
	l <- list()
	cl <- makeCluster(nc)
	registerDoParallel(cl)
	# read chunks of size n for each worker until EOF is reached
	while(length(x <- readLines(gfd, nc * n)) != 0){
		# if less is read, regenerate read indices to split by
		if(length(x) < nc * n)
			si <- splitIndices(length(x), nc)
		# split read for each worker
		x <- lapply(si, function(i) x[i])
		# apply transformation in parallel for each chunk
		ll <- foreach(x=x, sid=seq(sq, sq+nc-1), .inorder=FALSE, .multicombine=TRUE) %dopar%
			fn(x, sid)
		# append to list of results and increment sequence number for next run
		l <- c(l, ll)
		sq <- sq + nc
	}
	stopCluster(cl)
	close(gfd)
	# reorder chunks by sequence number and return a final data.frame with everything
	bind_rows(l[sort.int(sapply(l, function(x) x$sq[1]), index.return=TRUE)$ix])
}

# helper function for counting length of longest common prefix in a string
# vector
# x: string concatenation of overlapping gene indexes for a given chromosome
# g: gene vector for a given chromosome
# will return a vector of overlapping genes with common prefix
# the current threshold is a common prefix of size >= 1
cntpref <- function(x, g){
	# apply for each combination of overlapping genes
	sapply(x, function(x){
		# get index vector into g
		x <- as.integer(unlist(strsplit(x, ",")))
		# can't have a common prefix if less than two overlapping genes
		if(length(x) < 2)
			return("")
		# vector of overlapping gene names
		g <- as.character(g[x])
		# for each gene, get maximal length of common prefix, and keep
		# those that do have one > 0
		x <- x[sapply(g, function(x){
			# exclude self, compare to the rest
			l <- ifelse(g != x, nchar(g), 0)
			# try each prefix length starting from 1
			# l is the vector of name length and is subtracted from on each iteration
			# the loop is terminated once there are no remaining
			# strings or all are different
			# returned is the length of the longest common prefix
			# between x (reference) and any of the others prior to
			# termination
			for(i in 1:max(l)){
				y <- g[l > 0]
				if(all(ifelse(length(y)==0, TRUE, substr(x, 1, i) != substr(y, 1, i))))
					return(i-1)
				l <- l - 1
			}
			i
		}) > 0]
		# concatenate remaining genes and return
		paste0(x, collapse=",")
	})
}

# read and apply some corrections to the consensus gene list
mkrefseq <- function(f){
	read.table("hg19/ncbiRefSeq.txt.gz") %>%
		# rename interesting columns
		select(chr=V3, start=V5, end=V6, strand=V4, gene=V13) %>%
		# merge genes with same name on the same strand and chromosome
		group_by(chr, gene, strand) %>%
		summarize(start=min(start), end=max(end)) %>%
		# sort by chromosome, start and end
		arrange(chr, start, end) %>%
		# for each chromosome, generate a list of overlapping gene
		# intervals, then a column with overlapping genes with a common
		# name prefix
		group_by(chr) %>%
		mutate(over=sapply(1:length(start), function(x){
			paste0(which(
				start[x] >= start & start[x] <= end |
				end[x] >= start & end[x] <= end |
				end >= start[x] & end <= end[x] |
				start >= start[x] & start <= end[x]),
				collapse=",")}),
			overg=cntpref(over, gene)) %>%
		# mark overlapping genes with common prefix to merge: define
		# new name using prefix and a "_comb" suffix
		group_by(chr, overg, strand) %>%
		mutate(gene=ifelse(overg!="", paste0(as.character(gene[1]), "_", "comb"), as.character(gene))) %>%
		# merge marked genes
		group_by(chr, gene, strand) %>%
		summarize(start=min(start), end=max(end)) %>%
		ungroup %>%
		# sort and rearrange columns as needed by bedtools, using an
		# empty score column
		arrange(chr, start, end, strand) %>%
		mutate(score=rep(0)) %>%
		select(chr, start, end, gene, score, strand) %>%
		write.gzip(f)
}

# helper transform function for gzreadlines()
# takes a chunk of the gc5 wig file and the sequence number
# the output is a data.frame later bound together with all others. to speed up
# processing, we use use vectors for each column and preallocate them. we
# calculate the maximum number of rows we could possibly generate, and slice
# the vectors if we got less before returning
gc5cnt <- function(x, sq){
	# pessimistic estimated maximal length of 100kb bin rows to generate
	m <- length(x) / (100000 / 10)
	# preallocate column vectors
	chr <- character(m)
	start <- integer(m)
	end <- integer(m)
	# total gc count
	cnt <- integer(m)
	# total number of nucleotides in bin (some may be skipped)
	nt <- integer(m)
	# row index
	n <- 0
	# split input vector of strings by both whitespace and = (for both
	# possible line formats)
	x <- strsplit(x, "\\s|=")
	# state machine compressing wig rows
	for(i in x){
		# a new chromosome appears: end previous row and start count
		# over
		if(i[1] == "variableStep"){
			n <- n + 1
			chr[n] <- i[3]
			start[n] <- 0
			end[n] <- 100000
			cnt[n] <- 0
			nt[n] <- 0
		# this is the first row, so chromosome is unknown; use sequence
		# number for placeholder (will be corrected later), get start
		# position from first column and start a new row
		}else if(n == 0){
			n <- n + 1
			k <- as.integer((as.integer(i[1]) - 1) / 100000) * 100000
			chr[n] <- paste0("chrsq", sq)
			start[n] <- k
			end[n] <- k + 100000
			cnt[n] <- as.integer(i[2]) / 20
			nt[n] <- 5
		# we passed this bin's end boundary, start a new row
		}else if(as.integer(i[1]) - 1 > end[n]){
			n <- n + 1
			chr[n] <- chr[n-1]
			start[n] <- start[n-1] + 100000
			end[n] <- end[n-1] + 100000
			cnt[n] <- as.integer(i[2]) / 20
			nt[n] <- 5
		# continue counting for current row
		}else{
			cnt[n] <- cnt[n] + as.integer(i[2]) / 20
			nt[n] <- nt[n] + 5
		}
	}
	# return a named list with all parameters to be coerced into a data.frame later
	list(chr=chr[1:n], start=start[1:n], end=end[1:n], n=cnt[1:n], nt=nt[1:n], sq=rep(sq, n))
}

# the gc5base is a wig file with %gc over 5bp. it is prohibitively large
# uncompressed and cannot be stored in memory in one piece. the standard file
# format we convert everything to is a bed file (3 or more columns).
# gzreadlines will read the file in chunks and in parallel, combining
# information for 100kb bins based on position. since chunks will not
# necessarily contain enough information to infer which chromosome they are
# from, we need to rebuild this information afterwards. finally, we need to
# generate an actual bed file, which will be the %gc over 100kb bins (i.e. no
# further counting is done on this parameter)
mkgc5 <- function(f){
	df <- gzreadlines("hg19/hg19.gc5Base.wigVarStep.gz", 10*100000, gc5cnt)
	# fix chromosome names: chunks are ordered, so all chunks with a
	# sequence name placeholder before another chromosome name are actually
	# part of the previous chromosome
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
	# filter unused chromosomes and merge rows with same start and end (for
	# chunks that cut in the middle of a bin)
	chrom <- read.table("prep/hg19.txt", header=TRUE)
	df <- df %>%
		mutate(chr=as.character(l[chr])) %>%
		filter(chr %in% chrom$chrom) %>%
		group_by(chr, start, end) %>%
		summarize(n=sum(n), nt=sum(nt)) %>%
		ungroup
	# add more bins for missing rows at the end of each chromosome,
	# excluded from gc5 data
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
	# bind new rows with the rest, calculate actual %gc, create chr%02d
	# formatted chromosome names to guide sorting (like sort -V) and export
	df %>%
		rbind(df2) %>%
		mutate(gc=n/nt, chr0=ifelse(nchar(chr) == 4 & substr(chr, 4, 4) %in% 0:9, paste0("chr0", substr(chr, 4, 4)), chr)) %>%
		arrange(chr0, start) %>%
		select(chr, start, end, gc, n, nt, -chr0) %>%
		write.gzip(f)
}

# get hg19 chromosome size from ucsc mysql server, filtering uninteresting
# chromosomes, and sorting by name like sort -V
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

# extract a/b profile for each 100kb bin from excel table
mkab <- function(f){
	read_xlsx("gf/Table-AouBouAlways.xlsx", col_types="text") %>%
		mutate(start=format(as.integer(start)-1, scientific=FALSE, trim=TRUE)) %>%
		group_by(chr) %>%
		mutate(sign=sign(as.numeric(HUVEC)), sign=ifelse(sign==0,1,sign),
			trans=sign!=lead(sign), trans=ifelse(is.na(trans), FALSE, trans),
			flank=trans | lead(trans, default=FALSE) | lag(trans, default=FALSE) | lag(trans, 2, default=FALSE)) %>%
		ungroup %>%
		mutate(HUVECnoflank=ifelse(flank, NA, HUVEC)) %>%
		select(chr, start, end, AorBvec, HUVEC, HUVECnoflank, IMR90) %>%
		write.table(f, sep="\t", quote=FALSE, row.names=FALSE)
}

# extract repseq list and correlation from excel file
mkhuvecrep <- function(f){
	read_xlsx("gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx", col_types="text",
		.name_repair=~gsub(" ", "_", .x)) %>%
		select(-Rank_CorrB) %>%
		write.table(f, sep="\t", quote=FALSE, row.names=FALSE)
}

# concatenate groseq replicates for each strand
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

# extract silencer list from suppl table
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

# list all transformations to make with output filename
mkconv <- function(){
	l <- list(
		list(f="prep/hg19.refseq.bed.gz", fn=mkrefseq),
		list(f="prep/hg19.txt", fn=mkchr),
		list(f="prep/hg19.gc.bed.gz", fn=mkgc5),
		list(f="prep/ab.bed", fn=mkab),
		list(f="prep/huvec.repseq.tsv", fn=mkhuvecrep),
		list(f="prep/huvec.groseq.allrep.bed.gz", fn=mkgroseq),
		list(f="prep/huvec.silencer.bed.gz", fn=mksilencer)
	)
	for(i in l)
		if(file.access(i$f, 4) != 0)
			i$fn(i$f)
}

# split chromhmm track into beds for each class
mkchromhmm <- function(){
	l <- c("1_Active_Promoter",
		"2_Weak_Promoter",
		"3_Poised_Promoter",
		"4_Strong_Enhancer",
		"5_Strong_Enhancer",
		"6_Weak_Enhancer",
		"7_Weak_Enhancer")
	f <- paste0("prep/huvec.chromhmm.", tolower(sub("_", ".", l)), ".bed.gz")
	i <- which(file.access(f, 4) != 0)
	if(length(i) != 0){
		x <- read.table("hg19/wgEncodeBroadHmmHuvecHMM.bed.gz")
		x <- sapply(i, function(i){
			x %>%
				filter(V4 == l[i]) %>%
				mutate(V4=rep(sub("_.*", "", l[i]))) %>%
				select(V1, V2, V3, V4) %>%
				write.gzip(f[i])
			NULL
		})
	}
}

# generate beds for each subset of repseqs after extracting all repseqs from
# repeatmasker's output
mkrepseq <- function(){
	repmask <- read.table("hg19/hg19.fa.out.gz", skip=2) %>%
		select(chr=V5, start=V6, end=V7, id=V10)
	repcor <- read.table("prep/huvec.repseq.tsv", header=TRUE)
	psil <- read_xlsx("gf/Liste.ProtoSil.Forts-Etudiants-14nov2019.xlsx", col_types="text",
		.name_repair=~gsub(" ", "_", .x)) %>%
		select(4, 5) %>%
		mutate(prob=liste_FINALE, subb=liste) %>%
		filter(!row_number() %in% c(1,2))
	psilprob <- data.frame(Name=psil$prob)
	psilsubb <- data.frame(Name=psil$subb)

	l <- list(
		list(f="prep/huvec.repseq.te.bed.gz", q=quote(repcor %>% filter(Class %in% c("SINE", "RC", "SVA", "LTR", "LINE", "DNA")))),
		list(f="prep/huvec.repseq.ltr.bed.gz", q=quote(repcor %>% filter(Class == "LTR"))),
		list(f="prep/huvec.repseq.l1.bed.gz", q=quote(repcor %>% filter(Class == "LINE" & Family == "L1"))),
		list(f="prep/huvec.repseq.l1.ltr.bed.gz", q=quote(repcor %>% filter(Class %in% c("LTR", "LINE")))),
		list(f="prep/huvec.repseq.prob.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01))),
		list(f="prep/huvec.repseq.prob.te.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class %in% c("SINE", "RC", "SVA", "LTR", "LINE", "DNA")))),
		list(f="prep/huvec.repseq.prob.nonte.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class %in% c("tRNA", "snRNA", "Simple_repeat", "scRNA", "Satellite", "rRNA", "Low_complexity")))),
		list(f="prep/huvec.repseq.prob.sr.lc.sat.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class %in% c("Simple_repeat", "Low_complexity", "Satellite")))),
		list(f="prep/huvec.repseq.prob.ltr.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class == "LTR"))),
		list(f="prep/huvec.repseq.prob.l1.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class == "LINE" & Family == "L1"))),
		list(f="prep/huvec.repseq.prob.l1.ltr.bed.gz", q=quote(repcor %>% filter(Moyenne_corA < -0.01 & Class %in% c("LTR", "LINE")))),
		list(f="prep/huvec.repseq.psil.prob.bed.gz", q=quote(psilprob)),
		list(f="prep/huvec.repseq.proa.bed.gz", q=quote(repcor %>% filter(Moyenne_corA >= 0.01))),
		list(f="prep/huvec.repseq.proa.te.bed.gz", q=quote(repcor %>% filter(Moyenne_corA >= 0.01 & Class %in% c("SINE", "RC", "SVA", "LTR", "LINE", "DNA")))),
		list(f="prep/huvec.repseq.proa.nonte.bed.gz", q=quote(repcor %>% filter(Moyenne_corA >= 0.01 & Class %in% c("tRNA", "snRNA", "Simple_repeat", "scRNA", "Satellite", "rRNA", "Low_complexity")))),
		list(f="prep/huvec.repseq.trna.bed.gz", q=quote(repcor %>% filter(Class == "tRNA")))
	)
	f <- paste0("prep/repseq/huvec.repseq.only.",
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
	l <- foreach(f=f, n=n, .inorder=FALSE, .multicombine=TRUE, .export="write.gzip") %dopar%
		write.gzip(repmask[repmask$id %in% n,], f)
	stopCluster(cl)

	if(file.access("prep/huvec.repseq.l1long.bed.gz", 4) != 0)
		read.table("prep/huvec.repseq.l1.bed.gz") %>%
			filter(V3-V2 > 5000) %>%
			write.gzip("prep/huvec.repseq.l1long.bed.gz")
	if(file.access("prep/huvec.repseq.l1long.ltr.bed.gz", 4) != 0)
		read.table("prep/huvec.repseq.l1.ltr.bed.gz") %>%
			filter(V3-V2 > 5000) %>%
			write.gzip("prep/huvec.repseq.l1long.ltr.bed.gz")
	if(file.access("prep/huvec.repseq.prob.l1long.bed.gz", 4) != 0)
		read.table("prep/huvec.repseq.prob.l1.bed.gz") %>%
			filter(V3-V2 > 5000) %>%
			write.gzip("prep/huvec.repseq.prob.l1long.bed.gz")
	if(file.access("prep/huvec.repseq.prob.l1long.ltr.bed.gz", 4) != 0)
		read.table("prep/huvec.repseq.prob.l1.ltr.bed.gz") %>%
			filter(V3-V2 > 5000) %>%
			write.gzip("prep/huvec.repseq.prob.l1long.ltr.bed.gz")
}

# run all conversions
mkconv()
mkchromhmm()
mkrepseq()
