# generate linear models from counted parameters, diagnostic plots and summaries for each
library(dplyr)
library(broom)
library(ggplot2)
library(doParallel)
source("lib.R")

# generate model and model plots; takes the model data.frame, output directory
# and model name
# model is created on eigenvector with flanking regions filtered out, and is
# then run against raw eigenvector; bedgraphs for both are made
model <- function(ab, dir, name){
	# create model directory
	dir.create(dir)

	# list parameter names
	expl <- colnames(ab)[-c(1:5)]
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

	# get model diagnostic plots
	pdf(paste0(dir, "plots.pdf"), width=12.1, height=9.7)
	layout(matrix(c(1,2,3,4),2,2))
	plot(m)
	dev.off()

	# attempt for a better residuals vs fitted plot
	pdf(paste0(dir, "res.vs.fitted.pdf"), width=12.1, height=9.7)
	g <- data.frame(Fitted=m$fitted.values, Residuals=m$residuals) %>%
		ggplot(aes(Fitted, Residuals)) +
			geom_point(na.rm=TRUE, color="darkblue", alpha=0.2, shape=20) +
			geom_smooth(se=FALSE, color="red", size=0.4) +
			ggtitle("Model residuals versus fitted values")
	print(g)
	dev.off()

	# attempt for a better qqplot
	pdf(paste0(dir, "qqnorm.pdf"), width=12.1, height=9.7)
	g <- data.frame(Residuals=m$residuals) %>%
		ggplot(aes(sample=Residuals)) +
			stat_qq_line() +
			stat_qq(color="darkblue", size=0.7) +
			ggtitle("Quantile-quantile plot of residuals") +
			xlab("Theoretical normal distribution quantiles") +
			ylab("Standardized residuals")
	print(g)
	dev.off()

	# predict A/B profile on raw eigenvector; must use same variable names
	# if using predict()
	pred <- ab %>%
		mutate(eigenvectornf=eigenvector) %>%
		select(-eigenvector)
	v <- predict(m, pred)

	# attempt at a observed vs predicted plot
	pdf(paste0(dir, "obs.vs.fitted.pdf"), width=12.1, height=9.7)
	g <- pred %>%
		rename(Observed=eigenvectornf) %>%
		mutate(Fitted=v) %>%
		ggplot(aes(Fitted, Observed)) +
			geom_point(na.rm=TRUE, color="darkblue", alpha=0.2, shape=20) +
			geom_smooth(method="lm", na.rm=TRUE, se=FALSE, color="red", size=0.4) +
			ggtitle(paste0("Observed eigenvector values versus values predicted by model ",
				name, " (R²adj=", round(summary(m)$adj.r.squared, 2), ")"))
	print(g)
	dev.off()
}

# read table of all counts
ab <- read.table("tabs/lmparms.tsv.gz", header=TRUE)

# epigenomic and genomic parameter lists to cross-combine
# get list of all possible combinations between the two lists, and filenames of
# these combinations
l <- lapply(c("lm.eparm.tsv", "lm.gparm.tsv"), read.parms)
l <- comb.parms(l)
f <- paste0("lm/", names(l), "/")

# only generate missing models
i <- which(file.access(f, 4) != 0)
l <- l[i]
f <- f[i]

# generate a list of table slices for each combination, which foreach will
# distribute among workers
abl <- lapply(l, function(x){
	select(ab, eigenvector, eigenvectornf, !!!syms(as.character(x)))
})

# generate missing models
cl <- makeCluster(detectCores())
registerDoParallel(cl)
l <- foreach(ab=abl, f=f, n=names(l), .multicombine=TRUE, .inorder=FALSE, .packages=c("ggplot2", "dplyr", "broom")) %dopar% {
	model(ab, f, n)
	NULL
}
stopCluster(cl)

# read summaries of all models
l <- list.dirs("lm", full.names=FALSE)[-1]
perf <- lapply(l, function(x) read.table(paste0("lm/", x, "/perf.tsv"), header=TRUE))
coef <- lapply(l, function(x) cbind(read.table(paste0("lm/", x, "/coef.tsv"), header=TRUE), list(model=x)))

# write a summary of each model's performance, sorted by adjusted rsquare
bind_rows(perf) %>%
	mutate(model=l,
		parms=sapply(coef, function(x) paste(as.character(x[-1,1]), collapse=", "))) %>%
	select(model, everything()) %>%
	arrange(desc(adj.r.squared)) %>%
	write.tsv("lm/summary.tsv", col.names=TRUE)

# empty template data.frame for a table of parameters for each model
abz <- ab %>%
	select(-eigenvector, -eigenvectornf) %>%
	filter(rep(FALSE)) %>%
	mutate(model=character(), b0=numeric()) %>%
	select(model, b0, everything())
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
lapply(coeft, function(x) full_join(abz, x, by=colnames(x))) %>%
	bind_rows %>%
	write.tsv("lm/parms.tsv", col.names=TRUE)