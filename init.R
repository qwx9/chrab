# fetch input data from Encode, GEO, and others
library(dplyr)

# read a csv file describing each piece of data and extract file url and output
# file path
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

# check for required files that cannot be fetched
req <- c(
	file.path("gf", "Kassiotis-List.ORI.RepSeq.CorrB-Fourel.11july.xlsx"),
	file.path("gf", "Liste.ProtoSil.Forts-Etudiants-14nov2019.xlsx"),
	file.path("gf", "Table-AouBouAlways.xlsx")
)
for(i in req)
	if(file.access(i, 4) != 0)
		stop(paste("missing required file:", i))

# create all project directories if missing
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

# get data pointed to in each csv file
getfromcsv("encode.csv")
getfromcsv("geo.csv")
getfromcsv("hg19.csv")
getfromcsv("supp.csv")
