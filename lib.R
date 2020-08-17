# shared functions

# wrapper for write.table to output a tsv file
write.tsv <- function(x, con, col.names=FALSE, row.names=FALSE, ...){
	write.table(x, con, col.names=col.names, row.names=row.names, quote=FALSE, sep="\t", ...)
}

# wrapper for write.tsv to output a gz compressed file
write.gzip <- function(x, file, col.names=FALSE, ...){
	gfd <- gzfile(file, "wb")
	write.tsv(x, gfd, col.names=col.names, ...)
	close(gfd)
}

# write compressed bedgraph from data.frame with standard header
write.gzbedg <- function(x, file, name){
	gfd <- gzfile(file, "wb")
	cat(paste0("track type=bedGraph visibility=full",
		"color=200,100,0 altColor=0,100,200 priority=20 height=200 name=", name, "\n"),
		file=gfd)
	write.tsv(x, gfd, append=TRUE)
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
	if(!is.list(x) || length(x) < 2)
		stop("read.parms: must have at least two lists to combine")
	l <- apply(expand.grid(x[[1]], x[[2]]), 1, unlist)
	names(l) <- apply(expand.grid(seq_along(x[[1]])-1, seq_along(x[[2]])-1), 1, function(x){
		paste0("m", x[1], ".", x[2])
	})
	l[-1]
}

# slice existing files off a file path vector and associated vectors of same
# length (other metadata)
prune.extant <- function(df){
	df[which(file.access(df$f, 4) != 0),]
}
