write.gzip <- function(x, file, cn=FALSE, ...){
	gfd <- gzfile(file, "wb")
	write.table(x, gfd, col.names=cn, row.names=FALSE, quote=FALSE, sep="\t", ...)
	close(gfd)
}
