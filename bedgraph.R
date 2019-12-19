# take all generated counts and convert to bedgraphs
library(doParallel)

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
	# add element name to bedgraph header
	s <- sub("\\.bed\\.gz", "", sub("^(cnt/|cnt/repseq/)", "name=", l))
	# write a gz compressed file directly
	gfd <- gzfile(f, "wb")
	# write final bedgraph header
	cat(paste("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200", s, "\n"),
		file=gfd)
	# append count bed file contents
	write.table(read.table(l), gfd,
		col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	close(gfd)
}
stopCluster(cl)
