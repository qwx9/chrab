# install missing required packages, without reinstalling anything
l <- c("readxl", "dplyr", "ggplot2", "hexbin", "doParallel", "dbplyr", "RMySQL", "ggridges", "gridExtra", "tidyr", "broom")
l <- l[! l %in% installed.packages()[,"Package"]]
if(length(l) > 0)
	install.packages(l)
