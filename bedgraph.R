# take all generated counts and convert to bedgraphs
library(dplyr)
library(doParallel)
source("lib.R")

pullcnt <- function(ab, x){
	x <- read.table(x)
	if(ncol(x) > 1){
		x %>%
			filter(V1 %in% ab$chr) %>%
			select(4) %>%
			pull
	}else
		pull(x)
}

# get bin positions for files that don't contain them
ab <- read.table("prep/ab.bed", header=TRUE) %>%
	select(chr, start, end)

# read neural network results if they exist
if(file.access("gf/results.csv", 4) == 0){
	x <- read.table("gf/results.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	ab %>%
		mutate(pred=x$prediction_nn) %>%
		write.gzbedg("igv/nnpred.bedgraph.gz", "nnpred")
}

# generate list of files to generate along with vector of source tables and
# bedgraph names
# add count files
l <- list.files("cnt", pattern="*.gz", full.names=TRUE, recursive=TRUE)
l1 <- data.frame(l=l,
	f=sub("\\.bed", ".bedgraph", sub("^cnt/", "igv/", l)),
	n=gsub("^(cnt/|cnt/repseq/)|\\.bed\\.gz", "", l),
	stringsAsFactors=FALSE
)
# add model files
l <- list.files("score", pattern="pred.tsv.gz", full.names=TRUE, recursive=TRUE)
l2 <- data.frame(l=l,
	f=sub("\\.tsv", ".bedgraph", l),
	n=gsub("^score/|/pred.tsv.gz$", "", l),
	stringsAsFactors=FALSE
)

# don't overwrite existing files
l <- prune.extant(bind_rows(l1, l2))
if(nrow(l) == 0)
	quit()

nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
l <- foreach(l=l$l, f=l$f, n=l$n, .inorder=FALSE, .multicombine=TRUE, .packages="dplyr") %dopar% {
	ab %>%
		mutate(v=pullcnt(., l)) %>%
		filter(!is.na(v)) %>%
		write.gzbedg(f, n)
}
stopCluster(cl)
