# concatenate all counts into a single table, and generate other summary tables
suppressPackageStartupMessages({
library(dplyr)
})
source("lib.R")

# return a row of descriptive statistics for a variable in the filtered
# data.frame ab
statsrow <- function(ab, var){
	v <- ab[,var]
	su <- sum(v, na.rm=TRUE)
	if(all(is.na(v))){
		data.frame(parm=var,
			EanyA=NA,
			EanyB=NA,
			EalwaysA=NA,
			EalwaysB=NA,
			min=NA,
			max=NA,
			mean=NA,
			sd=NA,
			median=NA,
			q1=NA,
			q3=NA,
			na=length(v),
			stringsAsFactors=FALSE
		)
	}else{
		data.frame(parm=var,
			EanyA=sum(v[ab$anyA], na.rm=TRUE) / su * 100,
			EanyB=sum(v[ab$anyB], na.rm=TRUE) / su * 100,
			EalwaysA=sum(v[ab$alwaysA], na.rm=TRUE) / su * 100,
			EalwaysB=sum(v[ab$alwaysB], na.rm=TRUE) / su * 100,
			min=min(v, na.rm=TRUE),
			max=max(v, na.rm=TRUE),
			mean=mean(v, na.rm=TRUE),
			sd=sd(v, na.rm=TRUE),
			median=median(v, na.rm=TRUE),
			q1=quantile(v, 0.25, na.rm=TRUE),
			q3=quantile(v, 0.75, na.rm=TRUE),
			na=sum(is.na(v)),
			stringsAsFactors=FALSE
		)
	}
}

# cap a count based on tukey's fences
capcnt <- function(x, k=5){
	ubound <- quantile(x, 0.75) + k * IQR(x)
	ifelse(x > ubound, ubound, x)
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
write.counts <- function(x, file){
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
			write.tsv(paste0("tabs/", file, ".by", length(i), ".tsv"))
	})
}

# helper function to read a count table, filter data on unused chromosomes and
# pull count column, to add to a table
addcol <- function(chr, f){
	read.table(f) %>%
		filter(V1 %in% unique(chr)) %>%
		pull(V4)
}

# read initial binning table with eigenvector
ab <- read.table("prep/ab.bed", header=TRUE)
# list all count files (excluding those for individual repseqs), then add
# column for each, shortening the column name
for(i in list.files("cnt", pattern="*.gz", full.names=TRUE)){
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, i))
}

# special case: cap groseq score outliers to an upper boudary
ab <- ab %>%
	mutate(huvec.groseq.capped.meanofmean=capcnt(huvec.groseq.meanofmean),
		huvec.groseq.capped.meanofsum=capcnt(huvec.groseq.meanofsum),
		huvec.groseq.capped.sumofsum=capcnt(huvec.groseq.sumofsum))

# write out non-individual repseq count table
write.gzip(ab, "tabs/counts.tsv.gz", TRUE)

# epigenomic and genomic parameter lists for modeling
l <- lapply(c("lm.eparm.tsv", "lm.gparm.tsv"), read.parms)
# select modeling parameters and apply transformations
abm <- ab %>%
	select(HUVEC, HUVECnoflank, !!!syms(unique(unlist(l)))) %>%
	# normalize parameter columns range to [0;1] for easier interpretation
	# (doesn't affect prediction efficiency)
	mutate_if(1:ncol(.) > 2, ~. / max(., na.rm=TRUE)) %>%
	rename(eigenvector=HUVEC, eigenvectornf=HUVECnoflank)
# write out parameter values for neural network modeling
abm %>%
	select(-eigenvectornf) %>%
	write.csv(file="tabs/nnparms.csv", quote=FALSE, row.names=FALSE)
# write out parameter values for linear modeling
abm %>%
	write.gzip("tabs/lmparms.tsv.gz", TRUE)
# generate a correlation matrix for plotting
abm %>%
	select(-eigenvectornf) %>%
	cormat %>%
	write.tsv("tabs/lmparmscor.tsv", col.names=TRUE, row.names=TRUE)

# subdivide A/B regions into classes by gene density and presence of
# transcriptional activity, for each eigenvector (full and without flanking
# regions)
ab <- ab %>%
	mutate(class1=ifelse(HUVEC < 0, "B", "A"),
		class2=ifelse(hg19.refseq >= 4, "highgenedensity", ifelse(hg19.refseq > 0, "normalgenedensity", "nogene")),
		class2=factor(class2, levels=c("highgenedensity", "normalgenedensity", "nogene")),
		class3=ifelse(huvec.chromhmm.any.active.promoters > 0, "hasactive", "noactive"),
		class4=ifelse(HUVECnoflank < 0, "B", "A"),
		class=paste(class1, class2, AorBvec, class3, sep="_"),
		classF=paste(class4, class2, AorBvec, class3, sep="_"))
# export tables of counts split by each class
ab %>%
	select(type=class1, genedensity=class2, pc1=AorBvec, active=class3) %>%
	write.counts("bins")
ab %>%
	select(type=class4, genedensity=class2, pc1=AorBvec, active=class3) %>%
	write.counts("bins.noflank")
# export table of bin classification
ab %>%
	select(-class1, -class2, -class3, -class4) %>%
	select(chr, start, end, HUVEC, HUVECnoflank, class, classF) %>%
	write.gzip("tabs/class.tsv.gz", TRUE)

# read in columns for individual repseqs as before
for(i in list.files("cnt/repseq", pattern="*.gz", full.names=TRUE)){
	s <- gsub("\\.bed\\.gz$", "", gsub("^cnt/repseq/", "", i))
	ab <- ab %>%
		mutate(!!s:=addcol(chr, i))
}
# export table with counts for all parameters
ab %>%
	# remove unused columns, then reorder remaining
	select(-AorBvec, -IMR90, -matches("class[1-4]")) %>%
	select(chr, start, end, HUVEC, HUVECnoflank, class, classF, everything()) %>%
	write.gzip("tabs/aball.tsv.gz", TRUE)

# get list of parameter columns for stats
vars <- ab %>%
	select(-chr, -start, -end, -starts_with("class"), -AorBvec) %>%
	colnames

# add filtering variables for enrichment calculations
ab <- ab %>%
	mutate(anyA=class1=="A", anyB=class1=="B",
		alwaysA=AorBvec=="AlwaysA", alwaysB=AorBvec=="AlwaysB")

# list of subclass tables to generate with filter condition
l <- list(
	list(f="", q=quote(rep(TRUE))),
	list(f=".alwaysA.hasactive", q=quote(ab$AorBvec=="AlwaysA" & ab$huvec.chromhmm.any.active.promoters > 0)),
	list(f=".alwaysA.noactive", q=quote(ab$AorBvec=="AlwaysA" & ab$huvec.chromhmm.any.active.promoters == 0)),
	list(f=".alwaysB.hasactive", q=quote(ab$AorBvec=="AlwaysB" & ab$huvec.chromhmm.any.active.promoters > 0)),
	list(f=".alwaysB.noactive", q=quote(ab$AorBvec=="AlwaysB" & ab$huvec.chromhmm.any.active.promoters == 0)),
	list(f=".aorb.hasactive", q=quote(ab$AorBvec=="AorB" & ab$huvec.chromhmm.any.active.promoters > 0)),
	list(f=".aorb.noactive", q=quote(ab$AorBvec=="AorB" & ab$huvec.chromhmm.any.active.promoters == 0)),
	list(f=".alwaysA.highgene", q=quote(ab$AorBvec=="AlwaysA" & ab$class2 == "highgenedensity")),
	list(f=".alwaysA.reggene", q=quote(ab$AorBvec=="AlwaysA" & ab$class2 == "normalgenedensity")),
	list(f=".alwaysA.lowgene", q=quote(ab$AorBvec=="AlwaysA" & ab$class2 == "nogene")),
	list(f=".alwaysB.highgene", q=quote(ab$AorBvec=="AlwaysB" & ab$class2 == "highgenedensity")),
	list(f=".alwaysB.reggene", q=quote(ab$AorBvec=="AlwaysB" & ab$class2 == "normalgenedensity")),
	list(f=".alwaysB.lowgene", q=quote(ab$AorBvec=="AlwaysB" & ab$class2 == "nogene")),
	list(f=".aorba.highgene", q=quote(ab$AorBvec=="AorB" & ab$class1 == "A" & ab$class2 == "highgenedensity")),
	list(f=".aorba.reggene", q=quote(ab$AorBvec=="AorB" & ab$class1 == "A" & ab$class2 == "normalgenedensity")),
	list(f=".aorba.lowgene", q=quote(ab$AorBvec=="AorB" & ab$class1 == "A" & ab$class2 == "nogene")),
	list(f=".aorbb.highgene", q=quote(ab$AorBvec=="AorB" & ab$class1 == "B" & ab$class2 == "highgenedensity")),
	list(f=".aorbb.reggene", q=quote(ab$AorBvec=="AorB" & ab$class1 == "B" & ab$class2 == "normalgenedensity")),
	list(f=".aorbb.lowgene", q=quote(ab$AorBvec=="AorB" & ab$class1 == "B" & ab$class2 == "nogene"))
)
l <- c(l, lapply(unique(ab$class), function(x) list(f=paste0(".", x), q=bquote(ab$class == .(x)))))
# generate the subclass tables
l <- lapply(l, function(x){
	d <- filter(ab, eval(x$q))
	lapply(vars, function(n) statsrow(d, n)) %>%
		bind_rows %>%
		write.gzip(paste0("tabs/parms", x$f, ".tsv.gz"), TRUE)
})

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
