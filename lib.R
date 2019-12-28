# shared functions

# wrapper for write.table to output a gz compressed file
write.gzip <- function(x, file, cn=FALSE, ...){
	gfd <- gzfile(file, "wb")
	write.table(x, gfd, col.names=cn, row.names=FALSE, quote=FALSE, sep="\t", ...)
	close(gfd)
}

# read a list of parameter sublists; all combinations between two such lists
# will be generated afterwards
# to add mutually exclusive combinations with elements of just one of the two
# lists, prepend with an empty vector
read.parms <- function(file){
	append(list(character(0)), strsplit(readLines(file), "\t"))
}

# read two parameter tables and return all combinations between them
# the very first combination should be NULL, NULL and is removed
comb.parms <- function(x){
	if(length(x) < 2)
		stop("read.parms: must have at least two lists to combine")
	l <- apply(expand.grid(x[[1]], x[[2]]), 1, unlist)
	names(l) <- apply(expand.grid(seq_along(x[[1]])-1, seq_along(x[[2]])-1), 1, function(x){
		paste0("m", x[1], ".", x[2])
	})
	l[-1]
}
