write.gzip <- function(x, file, ...){
	gfd <- gzfile(file, "wb")
	write.table(x, gfd, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t", ...)
	close(gfd)
}
