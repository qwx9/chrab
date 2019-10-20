require(readxl)
require(dplyr)
options(error=recover)

# read encode data
lst <- read.csv("encode.csv")
# remove existing files
idx <- vector()
for(i in 1:nrow(lst))
	if(file.exists(as.character(lst[i,"file.path"])))
		idx <- c(idx, i)
lst <- lst[-idx,]
# download remaining files
if(nrow(lst) > 0)
	apply(lst, 1, function(x){ download.file(x["file.url"], x["file.path"]) })

# convert A/B profile excel table to tsv
x <- read_xlsx("gf/Table-AouBouAlways.xlsx", col_types="text") %>%
	select(chr, start, end, HUVEC, IMR90, AorBvec)
write.table(x, file="gf/ab.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# convert RepBase correlation with A/B profile excel table for HUVEC to tsv
x <- read_xlsx("gf/Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx", col_types="text") %>%
	select(-"Rank CorrB")
write.table(x, file="gf/huvec.repseq.tsv", sep="\t", quote=FALSE, row.names=FALSE)
