# generate linear models from counted parameters, diagnostic plots and summaries for each
suppressPackageStartupMessages({
library(dplyr)
library(tidyr)
library(broom)
library(doParallel)
})
source("lib.R")

# generate model and model plots; takes the model data.frame, output directory
# and model name
# model is created on eigenvector with flanking regions filtered out, and is
# then run against raw eigenvector; bedgraphs for both are made
model <- function(ab, dir){
	# create model directory
	dir.create(dir)

	# list parameter names
	expl <- colnames(ab)[-c(1:2)]
	# generate linear model formula
	fx <- quote(paste0("eigenvectornf ~", paste0(expl, collapse="+")))
	# compute model
	m <- lm(eval(fx), ab)

	# write out lm summary
	sink(paste0(dir, "summary.txt"))
	print(summary(m))
	sink()

	# write tidy coefficient estimations, model performance statistics and
	# predictions
	tidy(m) %>%
		write.tsv(paste0(dir, "coef.tsv"), col.names=TRUE)
	glance(m) %>%
		write.tsv(paste0(dir, "perf.tsv"), col.names=TRUE)
	augment(m, newdata=ab) %>%
		select(.fitted) %>%
		write.gzip(paste0(dir, "pred.tsv.gz"))
}

# list all lm parameter tables
f <- list.files("tabs", pattern="\\.lm\\.tsv\\.gz", full.names=TRUE)
l <- lapply(f, function(x){
	# cell line name
	cell <- sub("tabs/", "", x) %>%
		sub(pattern="\\.lm\\.tsv\\.gz", replacement="")
	# read table of all counts
	ab <- read.table(x, header=TRUE)

	# epigenomic and genomic parameter lists to cross-combine
	# get list of all possible combinations between the two lists, and filenames of
	# these combinations
	l <- lapply(c(paste0(cell, ".lmparm.tsv"), "hg19.lmparm.tsv"), read.parms)
	l <- comb.parms(l)
	l$mall <- unique(unlist(l))
	f <- paste0("lm/", cell, "/", names(l), "/")
	dir.create(paste0("lm/", cell), showWarnings=FALSE)

	# only generate missing models
	i <- which(file.access(f, 4) != 0)
	l <- l[i]
	f <- f[i]
	if(length(i) == 0)
		return(NULL)

	# generate a list of table slices for each combination, which foreach will
	# distribute among workers
	abl <- lapply(l, function(x){
		select(ab, eigenvector, eigenvectornf, !!!syms(as.character(x)))
	})
	names(abl) <- f
	abl
}) %>%
	unlist(recursive=FALSE)
if(length(l) == 0)
	quit()

# generate missing models
cl <- makeCluster(detectCores())
registerDoParallel(cl)
l <- foreach(ab=l, f=names(l), .multicombine=TRUE, .inorder=FALSE, .packages=c("ggplot2", "dplyr", "broom")) %dopar% {
	model(ab, f)
	NULL
}
stopCluster(cl)

# get cell line names again
cell <- lapply(f, function(x) sub("tabs/", "", x) %>% sub(pattern="\\.lm\\.tsv\\.gz", replacement=""))
for(i in cell){
	pref <- paste0("lm/", i, "/")
	# read summaries of all models
	l <- list.dirs(pref, full.names=FALSE)[-1]
	perf <- lapply(l, function(x) read.table(paste0(pref, x, "/perf.tsv"), header=TRUE))
	coef <- lapply(l, function(x) cbind(read.table(paste0(pref, x, "/coef.tsv"), header=TRUE), list(model=x)))

	# write a summary of each model's performance, sorted by adjusted rsquare
	bind_rows(perf) %>%
		mutate(model=l,
			parms=sapply(coef, function(x) paste(as.character(x[-1,1]), collapse=", "))) %>%
		select(model, everything()) %>%
		arrange(desc(adj.r.squared)) %>%
		write.tsv(paste0(pref, "summary.tsv"), col.names=TRUE)

	# empty template data.frame for a table of parameters for each model
	# read parameter lists again
	l <- lapply(c(paste0(i, ".lmparm.tsv"), "hg19.lmparm.tsv"), read.parms)
	nm <- c("model", "b0", unique(unlist(l)))
	m <- matrix(rep(NA), nrow=0, ncol=length(nm), byrow=TRUE)
	colnames(m) <- nm
	m <- as.data.frame(m) %>%
		mutate_each(as.numeric) %>%
		mutate_at("model", as.character)

	# tidyverse "simple" transpose
	coeft <- lapply(coef, function(x){
		x %>%
			gather(var, val, 2:ncol(x), -model) %>%
			spread("term", "val", "model", convert=TRUE) %>%
			rename(b0="(Intercept)") %>%
			filter(var == "estimate") %>%
			select(-var) %>%
			mutate_if(is.factor, as.character)
	})
	# export a table for parameter used across all models
	lapply(coeft, function(x) full_join(m, x, by=colnames(x))) %>%
		bind_rows %>%
		cbind(bind_rows(perf) %>% select(adj.r.squared)) %>%
		select(model, adj.r.squared, everything()) %>%
		write.tsv(paste0(pref, "models.tsv"), col.names=TRUE)
}
