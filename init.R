require(readxl)
require(dplyr)
options(error=recover)

getfromcsv <- function(file){
	lst <- read.csv(file)

	# remove existing files
	idx <- vector()
	for(i in 1:nrow(lst))
		if(file.exists(as.character(lst[i,"file.path"])))
			idx <- c(idx, i)
	if(length(idx) > 0)
		lst <- lst[-idx,]

	# download remaining files
	if(nrow(lst) > 0)
		apply(lst, 1, function(x){ download.file(x["file.url"], x["file.path"]) })
}

getfromcsv("encode.csv")
getfromcsv("geo.csv")

# convert A/B profile excel table to tsv
x <- read_xlsx("gf/Table-AouBouAlways.xlsx", col_types="text") %>%
	select(chr, start, end, AorBvec, HUVEC, IMR90) %>%
	mutate(start=format(as.integer(start)-1, scientific=FALSE, trim=TRUE))
write.table(x, file="gf/ab.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# convert RepBase correlation with A/B profile excel table for HUVEC to tsv
x <- read_xlsx("gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx", col_types="text") %>%
	select(-"Rank CorrB")
write.table(x, file="gf/huvec.repseq.tsv", sep="\t", quote=FALSE, row.names=FALSE)
