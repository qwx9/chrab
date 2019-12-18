library(dplyr)

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

req <- c(
	file.path("gf", "Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx"),
	file.path("gf", "Table-AouBouAlways.xlsx")
)
for(i in req)
	if(file.access(i, 4) != 0)
		stop(paste("missing required file:", i))

dirs <- c(
	"hg19",
	"huvec",
	"imr90",
	"prep",
	"cnt",
	"plot",
	"tabs",
	"score",
	"prep/repseq",
	"cnt/repseq",
	"plot/repseq",
	"igv/repseq"
)
for(i in dirs)
	dir.create(i, showWarnings=FALSE)

getfromcsv("encode.csv")
getfromcsv("geo.csv")
getfromcsv("hg19.csv")
getfromcsv("supp.csv")
