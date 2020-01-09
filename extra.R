library(dplyr)
source("lib.R")

l <- list.dirs("lm", full.names=FALSE)[-1]
ab <- read.table("prep/ab.bed", header=TRUE) %>%
	select(chr, start, end)
d <- lapply(l, function(x){
	read.table(paste0("lm/", x, "/pred.bedgraph.gz"), skip=1, col.names=c("chr", "start", "end", x))
})
for(i in seq_along(d))
	ab <- ab <- merge(ab, d[[i]], by=c("chr", "start", "end"), all=TRUE)
write.gzip(ab, "extra/lm.predictions.tsv.gz", col.names=TRUE)
