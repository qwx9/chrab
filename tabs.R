# concatenate all counts into a single table, and generate other summary tables
suppressPackageStartupMessages({
library(dplyr)
})
source("lib.R")

enrichment <- function(v, class){
	vratio <- sum(v[class], na.rm=TRUE) / sum(v, na.rm=TRUE)
	binratio <- sum(class, na.rm=TRUE) / length(na.omit(class))
	vratio / binratio
}

# return a row of descriptive statistics for a variable in the filtered
# data.frame ab
statsrow <- function(ab, var, pc1list){
	v <- ab[,var]
	ispc1 <- var %in% pc1list
	if(all(is.na(v))){
		data.frame(parm=var,
			EanyA=NA,
			EanyB=NA,
			EalwaysA=NA,
			EalwaysB=NA,
			sum=NA,
			min=NA,
			max=NA,
			mean=NA,
			sd=NA,
			median=NA,
			q1=NA,
			q3=NA,
			na=length(v)
		)
	}else{
		# don't calculate enrichment and things that don't make sense for eigenvectors
		data.frame(parm=var,
			EanyA=ifelse(ispc1, NA, enrichment(v, ab$anyA)),
			EanyB=ifelse(ispc1, NA, enrichment(v, ab$anyB)),
			EalwaysA=ifelse(ispc1, NA, enrichment(v, ab$alwaysA)),
			EalwaysB=ifelse(ispc1, NA, enrichment(v, ab$alwaysB)),
			sum=ifelse(ispc1, NA, sum(v, na.rm=TRUE)),
			min=min(v, na.rm=TRUE),
			max=max(v, na.rm=TRUE),
			mean=mean(v, na.rm=TRUE),
			sd=sd(v, na.rm=TRUE),
			median=median(v, na.rm=TRUE),
			q1=quantile(v, 0.25, na.rm=TRUE),
			q3=quantile(v, 0.75, na.rm=TRUE),
			na=sum(is.na(v))
		)
	}
}

# set var-related columns to log10 scale
tolog10 <- function(tab, cell, var, cols){
	for(i in cols){
		new <- paste(cell, var, "log10", i, sep=".")
		old <- paste(cell, var, i, sep=".")
		tab <- tab %>%
			mutate(!!new:=log10(!!sym(old)+1)) %>%
				select(-!!sym(old))
	}
	tab
}

# generate a clustered correlation matrix
cormat <- function(x){
	# calculate correlation between all parameters, removing NAs
	x <- cor(x, use="complete.obs")
        # perform hierarchical clustering using the correlation as a
        # dissimilarity measure and order correlation matrix by clusters
        hc <- hclust(as.dist((1-x)/2))
        x[hc$order, hc$order]
}

# helper function for writing out counts split by one or more classes
write.counts <- function(x, file, cell){
	l <- colnames(x)
	# split iteratively over more and more subclasses and export table
	x <- sapply(seq_along(l), function(i){
		if(i > 1){
			i <- 1:i
			l <- l[i]
		}else{
			l <- "."
		}
		x %>%
			select(i) %>%
			# count NAs
			table(useNA="always") %>%
			as.data.frame %>%
			# rename automatically named Freq column
			rename(count=Freq) %>%
			# arrange as we specified
			arrange(!!!syms(l)) %>%
			write.tsv(paste0("tabs/", cell, ".", file, ".by", length(i), ".tsv"))
	})
}

# read a list of count tables, retain only element counts,
# bind all as columns of a dataframe
readtabs <- function(files, chr){
	lapply(files, function(x){
		read.table(x) %>%
		filter(V1 %in% chr) %>%
		select(!!gsub("\\.bed\\.gz$", "", gsub(".*/", "", x)):=V4)
	}) %>%
		bind_cols
}

mkcountsonly <- function(df, cell, abref){
	# read binning table from another cell line to select chr
	ab <- read.table(paste0("prep/", cell, ".ab.bed"), header=TRUE)
	# add cell counts and genomic elements
	l <- list.files("cnt", pattern=cell, full.names=TRUE)
	ab <- cbind(ab, readtabs(l, unique(ab$chr))) %>%
		cbind(df %>% filter(chr %in% unique(ab$chr)) %>% select(-chr)) %>%
		select(-AorBvec, -5:6) %>%
		# special case: log10 scale for h3k9me3
		tolog10(cell, "h3k9me3", c("mean", "sum"))
	# export count table with classes
	ab %>%
		write.gzip(paste0("tabs/", cell, ".countsbybin.tsv.gz"), TRUE)

	# get list of parameter columns for stats
	vars <- ab %>%
		select(-chr, -start, -end) %>%
		colnames
	# export global stats table without classes
	lapply(vars, function(n) statsrow(ab, n, n)) %>%
		bind_rows %>%
		write.gzip(paste("tabs/stats.global", cell, "tsv.gz", sep="."), TRUE)
}

# generate tables for each cell type
mktab <- function(df, cell, pc1, pc1nf){
	# read initial binning table with eigenvector
	ab <- read.table(paste0("prep/", cell, ".ab.bed"), header=TRUE)

	# add cell counts and genomic elements
	l <- list.files("cnt", pattern=cell, full.names=TRUE)
	ab <- cbind(ab, readtabs(l, unique(ab$chr))) %>%
		cbind(df %>% filter(chr %in% unique(ab$chr)) %>% select(-chr))

	# special case: log10 scale for huvec groseq
	if(cell == "huvec")
		ab <- tolog10(ab, cell, "groseq", c("meanofmean", "meanofsum", "sumofsum"))
	# special case: log10 scale for h3k9me3
	ab <- tolog10(ab, cell, "h3k9me3", c("mean", "sum"))

	# get list of parameter columns for stats
	vars <- ab %>%
		select(-chr, -start, -end, -AorBvec) %>%
		colnames

	# column indicating active promoters
	act <- paste0(cell, ".chromhmm.any.active.promoters")

	# subdivide A/B regions into classes by gene density and presence of
	# transcriptional activity, for each eigenvector (full and without flanking
	# regions)
	ab <- ab %>%
		mutate(class1=ifelse(!!sym(pc1) < 0, "B", "A"),
			class2=ifelse(hg19.refseq >= 4, "highgenedensity", ifelse(hg19.refseq > 0, "normalgenedensity", "nogene")),
			class2=factor(class2, levels=c("highgenedensity", "normalgenedensity", "nogene")),
			class3=ifelse(!!sym(act) > 0, "hasactive", "noactive"),
			class4=ifelse(!!sym(pc1nf) < 0, "B", "A"),
			class=paste(class1, class2, AorBvec, class3, sep="_"),
			classF=paste(class4, class2, AorBvec, class3, sep="_"))
	# export tables of counts split by each class
	ab %>%
		select(type=class1, genedensity=class2, pc1=AorBvec, active=class3) %>%
		write.counts("bincountsbyclass", cell)
	ab %>%
		select(type=class4, genedensity=class2, pc1=AorBvec, active=class3) %>%
		write.counts("bincountsbyclass.noflank", cell)
	# export table of bin classification
	ab %>%
		select(chr, start, end, !!sym(pc1), !!sym(pc1nf), class, classF) %>%
		write.gzip(paste0("tabs/", cell, ".classbybin.tsv.gz"), TRUE)

	# export count table with classes
	ab %>%
		# remove unused columns, then reorder remaining
		select(-matches("class[1-4]")) %>%
		select(chr, start, end, !!sym(pc1), class, !!sym(pc1nf), classF, everything()) %>%
		write.gzip(paste0("tabs/", cell, ".countsbybin.tsv.gz"), TRUE)

	# create binary count beds for each class generated above (for diagnostics)
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

	# generate distribution statistics tables
	# add filtering variables for enrichment calculations
	ab <- ab %>%
		mutate(anyA=class1=="A", anyB=class1=="B",
			alwaysA=AorBvec=="AlwaysA", alwaysB=AorBvec=="AlwaysB")
	# export global stats table without classes
	lapply(vars, function(n) statsrow(ab, n, c(pc1, pc1nf))) %>%
		bind_rows %>%
		write.gzip(paste("tabs/stats.global", cell, "tsv.gz", sep="."), TRUE)
	# stats tables for each individual class
	# list of subclass tables to generate with filter condition: these are selected classes constructed for later
	l <- list(
		list(f="alwaysA.hasactive", q=quote(AorBvec=="AlwaysA" & !!sym(act) > 0)),
		list(f="alwaysA.noactive", q=quote(AorBvec=="AlwaysA" & !!sym(act) == 0)),
		list(f="alwaysB.hasactive", q=quote(AorBvec=="AlwaysB" & !!sym(act) > 0)),
		list(f="alwaysB.noactive", q=quote(AorBvec=="AlwaysB" & !!sym(act) == 0)),
		list(f="aorb.hasactive", q=quote(AorBvec=="AorB" & !!sym(act) > 0)),
		list(f="aorb.noactive", q=quote(AorBvec=="AorB" & !!sym(act) == 0)),
		list(f="alwaysA.highgene", q=quote(AorBvec=="AlwaysA" & class2 == "highgenedensity")),
		list(f="alwaysA.reggene", q=quote(AorBvec=="AlwaysA" & class2 == "normalgenedensity")),
		list(f="alwaysA.lowgene", q=quote(AorBvec=="AlwaysA" & class2 == "nogene")),
		list(f="alwaysB.highgene", q=quote(AorBvec=="AlwaysB" & class2 == "highgenedensity")),
		list(f="alwaysB.reggene", q=quote(AorBvec=="AlwaysB" & class2 == "normalgenedensity")),
		list(f="alwaysB.lowgene", q=quote(AorBvec=="AlwaysB" & class2 == "nogene")),
		list(f="aorba.highgene", q=quote(AorBvec=="AorB" & class1 == "A" & class2 == "highgenedensity")),
		list(f="aorba.reggene", q=quote(AorBvec=="AorB" & class1 == "A" & class2 == "normalgenedensity")),
		list(f="aorba.lowgene", q=quote(AorBvec=="AorB" & class1 == "A" & class2 == "nogene")),
		list(f="aorbb.highgene", q=quote(AorBvec=="AorB" & class1 == "B" & class2 == "highgenedensity")),
		list(f="aorbb.reggene", q=quote(AorBvec=="AorB" & class1 == "B" & class2 == "normalgenedensity")),
		list(f="aorbb.lowgene", q=quote(AorBvec=="AorB" & class1 == "B" & class2 == "nogene"))
	)
	# add all the other class combinations as well
	l <- c(l, lapply(unique(ab$class), function(x) list(f=paste0("old.", x), q=bquote(ab$class == .(x)))))
	# generate the subclass tables
	l <- lapply(l, function(x){
		d <- filter(ab, eval(x$q))
		lapply(vars, function(n) statsrow(d, n, c(pc1, pc1nf))) %>%
			bind_rows %>%
			write.gzip(paste(paste0("tabs/stats/", cell), "statsbyclass", x$f, "tsv.gz", sep="."), TRUE)
	})

	# export adjusted counts for linear modeling
	# epigenomic and genomic parameter lists for modeling
	l <- lapply(paste(c(cell, "hg19"), "lmparm.tsv", sep="."), read.parms) %>%
		unlist %>%
		unique
	# select modeling parameters and apply transformations
	abm <- ab %>%
		select(!!sym(pc1), !!sym(pc1nf), !!!syms(l)) %>%
		# normalize parameter columns range to [0;1] for easier interpretation
		# (doesn't affect prediction efficiency)
		mutate_if(1:ncol(.) > 2, ~. / max(., na.rm=TRUE)) %>%
		rename(eigenvector=!!sym(pc1), eigenvectornf=!!sym(pc1nf))
	# write out parameter values for linear modeling
	abm %>%
		write.gzip(paste0("tabs/", cell, ".lm.tsv.gz"), TRUE)
	# generate a correlation matrix for plotting
	abm %>%
		select(-eigenvectornf) %>%
		cormat %>%
		write.tsv(paste0("tabs/", cell, ".lmcor.tsv"), col.names=TRUE, row.names=TRUE)

	# write out parameter values for neural network modeling
	abm %>%
		select(-eigenvectornf) %>%
		write.csv(file=paste0("tabs/", cell, ".nnet.csv"), quote=FALSE, row.names=FALSE)
}

dir.create("tabs/stats", showWarnings=FALSE)

# get chromosome names for each bin to filter later
chr <- read.table("prep/hg19w.bed") %>%
	select(chr=V1)
# read in genomic elements to reuse in each cell type tables
l <- c(
	list.files("cnt", pattern="hg19.*gz", full.names=TRUE),
	list.files("cnt/repseq", full.names=TRUE)
)
df <- cbind(chr, readtabs(l, unique(chr$chr)))
mktab(df, "huvec", "HUVEC", "HUVECnoflank")
mktab(df, "gm12878", "GM12878", "GM12878noflank")
mkcountsonly(df, "bcell")
