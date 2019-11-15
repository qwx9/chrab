l <- c("readxl", "dplyr", "ggplot2", "hexbin", "doParallel", "dbplyr", "RMySQL", "ggridges", "gridExtra")
l <- l[! l %in% installed.packages()[,"Package"]]
if(length(l) > 0)
	install.packages(l)
