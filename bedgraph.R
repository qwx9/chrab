require(doParallel)

l <- list.files("cnt", pattern="*.gz", full.names=TRUE, recursive=TRUE)
f <- sub("\\.bed", ".bedgraph", sub("^cnt/", "igv/", l))
i <- which(file.access(f, 4) != 0)
if(length(i) != 0){
	l <- l[i]
	f <- f[i]
}
nc <- detectCores()
cl <- makeCluster(nc)
registerDoParallel(cl)
l <- foreach(l=l, f=f, .inorder=FALSE) %dopar% {
	s <- sub("\\.bed\\.gz", "", sub("^(cnt/|cnt/repseq/)", "name=", l))
	gfd <- gzfile(f, "wb")
	cat(paste("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200", s, "\n"),
		file=gfd)
	write.table(read.table(l), gfd,
		col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	close(gfd)
}
stopCluster(cl)
