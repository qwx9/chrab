# take all generated counts and convert to bedgraphs
library(dplyr)
library(doParallel)
source("lib.R")

# get bin positions for files that don't contain them
ab <- read.table("prep/ab.bed", header=TRUE, stringsAsFactors=FALSE) %>%
	select(chr, start, end)

# read neural network results if they exist
if(file.access("gf/results.csv", 4) == 0){
	x <- read.table("gf/results.csv", sep="\t", header=TRUE, stringsAsFactors=FALSE)
	ab %>%
		mutate(pred=x$prediction_nn) %>%
		write.gzbedg("igv/nnpred.bedgraph.gz", "nnpred")
}

# recursively list count files, use same path for destination, and check if the
# bedgraphs already exist
l <- list.files("cnt", pattern="*.gz", full.names=TRUE, recursive=TRUE)
f <- sub("\\.bed", ".bedgraph", sub("^cnt/", "igv/", l))
i <- which(file.access(f, 4) != 0)
# no new files to be added
if(length(i) == 0)
	quit()
# generate only missing files, prep for parallelization by taking input data
# which foreach will distribute to workers on its own
l <- l[i]
f <- f[i]
nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
l <- foreach(l=l, f=f, .inorder=FALSE, .multicombine=TRUE) %dopar% {
	write.gzbedg(read.table(l), f, gsub("^(cnt/|cnt/repseq/)|\\.bed\\.gz", "", l))
}
stopCluster(cl)
