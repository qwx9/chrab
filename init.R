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
