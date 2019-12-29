# concatenate all counts into a single table, and generate other summary tables
library(dplyr)
source("lib.R")

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
# special case: cap groseq.score outliers to an upper boudary defined by
# tukey's fences with k=5, generating a new column
cap <- quantile(ab$huvec.groseq.score, 0.75) + 5 * IQR(ab$huvec.groseq.score)
ab <- ab %>%
        mutate(huvec.groseq.capped=ifelse(huvec.groseq.score > cap, cap, huvec.groseq.score))

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
# remove helper subclasses
ab <- ab %>%
	select(-class1, -class2, -class3, -class4)
# export table of bin classification
ab %>%
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
	select(-AorBvec, -IMR90) %>%
	select(chr, start, end, HUVEC, HUVECnoflank, class, classF, everything()) %>%
	write.gzip("tabs/aball.tsv.gz", TRUE)

# export descriptive statistics table for each count (removing non-count columns)
n <- colnames(ab)
n <- n[grep("^(NA_|[AB]_|AorB|class|chr|start|end)", n, invert=TRUE)]
# bind together rows for each parameter
bind_rows(lapply(n, function(n){
	x <- ab[,n]
	y <- na.omit(x)
	data.frame(min=min(y), max=max(y), mean=mean(y), sd=sd(y), median=median(y),
		q1=quantile(y, 0.25), q3=quantile(y, 0.75), na=sum(is.na(x)))
})) %>%
	# reorder columns to have parameter first
	mutate(parm=n) %>%
	select(parm, everything()) %>%
	write.gzip("tabs/parms.tsv.gz", TRUE)

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
