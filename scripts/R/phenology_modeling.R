# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2022.12.14

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Jotr_phenology-main")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("embarcadero")
library("ggdark")

source("../shared/Rscripts/base_graphics.R")

#-----------------------------------------------------------
# initial file loading

# flow <- read.csv("output/flowering_obs_climate.csv") # flowering/not flowering, gridded and annualized
# flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized and normed

flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)

ggplot(flow2, aes(x=lon, y=lat, color=flr)) + geom_point() + facet_wrap("year") + theme_bw()


# variant datasets -- dealing with the second flowering in 2019
flow2 <- flow %>% filter(!(year==2019.5 & flr==TRUE)) %>% mutate(year=floor(year)) # drop the late-flowering anomaly
flow3 <- flow
flow3$year[flow3$year==2019.5] <- 2019 # or merge 2019.5 into 2019?

glimpse(flow2) # 2,503 in our final working set

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="YUJA") 
yubr <- filter(flow2, type=="YUBR")

glimpse(yuja) # 1,096 obs (after iffy ones excluded)
glimpse(yubr) # 1,339 obs

#-------------------------------------------------------------------------
# fit candidate BART models, stepwise

# predictors
xnames <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1", "vpdmaxW0", "vpdminW0", "vpdmaxW0vW1", "vpdminW0vW1") # year + weather data, "curated"

# Full range ------------------------------------

# PCA to describe predictor variation ------

jotrPCA <- prcomp(flow2[,xnames], scale.=TRUE)
summary(jotrPCA)

# what do we get from this ...?
expc <- function(pca, ax){

	ord = order(abs(pca$rotation[,ax]),decreasing=TRUE)

	return(list(prop.var=round(100*summary(pca)$importance[2,ax]), rots=data.frame(Var=rownames(pca$rotation)[ord], Rot=pca$rotation[ord,ax])))
}
	
expc(jotrPCA, 1) # 42% of variance; higher values => lower pptY0
expc(jotrPCA, 2) # 24% of variance; higher values => greater vpdmaxW0vW1 (vpdmaxW0 > vpdmaxW1)
expc(jotrPCA, 3) # 14% of variance; higher values => lower tminW0vW1 (tminW0 < tminW1)
expc(jotrPCA, 4) # 9% of variance; higher values => greater tminW0

flow2PCA <- cbind(flow2,jotrPCA$x[,1:4])

ggplot(flow2PCA, aes(x=-PC1, y=PC2, color=flr)) + geom_point(alpha=0.5) + theme_bw()
ggplot(flow2PCA, aes(x=PC2, y=PC3, color=flr)) + geom_point(alpha=0.5) + theme_bw()
ggplot(flow2PCA, aes(x=PC3, y=-PC1, color=flr)) + geom_point(alpha=0.5) + theme_bw()
ggplot(flow2PCA, aes(x=-PC1, y=PC4, color=flr)) + geom_point(alpha=0.5) + theme_bw()



jotrPCAvarimp <- varimp.diag(y.data=as.numeric(flow2PCA$flr), x.data=flow2PCA[,c("PC1","PC2","PC3","PC4")], ri.data=flow2PCA$year)


# variable importance across the whole predictor set

# doing replicate subsamples to cope with unbalanced data ...
# BY YEAR ...
table(flow2$year, flow2$flr) # so that's fun; the imbalance is not consistent
flow2ssNs <- apply(table(flow2$year, flow2$flr), 1, min) # what's the smaller sample size in each year

nreps <- 10 # start small

flow2subs <- vector("list", nreps)

# okay let's try this ...
for(ss in 1:nreps){ # loop over subsample replicates

subsam <- NULL

	for(yr in unique(flow2$year)){ # loop over years to balance that way ...
	
	subN <- min(table(flow2[flow2$year==yr,]$flr))
	subsam <- flow2 %>% filter(year==yr) %>% group_by(flr) %>% slice_sample(n=subN) %>% rbind(subsam, .)
	
	}

flow2subs[[ss]] <- subsam

}

class(flow2subs)
head(flow2subs[[1]])

save(flow2subs, file="output/BART/jotr.flowering.obs_subsamples-by-year.Rdata") # save this for replicability and eventually model-fitting ...


jotr.varimps <- lapply(flow2subs, function(ss) varimp.diag(y.data=as.numeric(ss$flr), x.data=ss[,xnames], ri.data=ss$year)) # fingers crossed

{cairo_pdf(file = "output/figures/jotr.subsamp-by-year.varimps.pdf", width=6, height=4, onefile=TRUE)

jotr.varimps

}
dev.off() # OKAY this worked as intended

# ACROSS YEARS ...
min(table(flow2$flr)) # what's the smaller sample size ah right flowering

nreps <- 10 # start small

flow2subs <- vector("list", nreps)

# okay let's try this ...
for(ss in 1:nreps){ # loop over subsample replicates

flow2subs[[ss]] <- flow2 %>% group_by(flr) %>% slice_sample(n=min(table(flow2$flr)))

}

class(flow2subs)
head(flow2subs[[1]])

save(flow2subs, file="output/BART/jotr.flowering.obs_subsamples-across-years.Rdata") # save this for replicability and eventually model-fitting ...


jotr.varimps <- lapply(flow2subs, function(ss) varimp.diag(y.data=as.numeric(ss$flr), x.data=ss[,xnames], ri.data=ss$year)) # fingers crossed

{cairo_pdf(file = "output/figures/jotr.subsamp-across-year.varimps.pdf", width=6, height=4, onefile=TRUE)

jotr.varimps

}
dev.off() # OKAY this worked as intended


# stepwise model fitting
jotr.flr.mod.step <- bart.step(y.data=as.numeric(flow2[,"flr"]), x.data=flow2[,xnames], ri.data=flow2[,"year"], full=FALSE, quiet=TRUE) 
save(jotr.flr.mod.step, file="output/BART/bart.step.models.jotr.Rdata")

# load(file="output/BART/bart.step.models.jotr.Rdata")

summary(jotr.flr.mod.step)

varimp(jotr.flr.mod.step, plots=TRUE)
partial(jotr.flr.mod.step) # this needs debugging, per CJC
plot(jotr.flr.mod.step) #??

# refitting with vars indicated by varimp:
jotr.preds <- c("tmaxW0vW1", "pptY0", "tmaxW0", "tminW0")

jotr.mod <- rbart_vi(as.formula(paste(paste('flr', paste(jotr.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = flow2,
	group.by = flow2[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(jotr.mod)
varimp(jotr.mod)

save(jotr.mod, file="output/BART/bart.ri.model.Jotr.Rdata")

# load(file="output/BART/bart.ri.model.Jotr.Rdata")

summary(jotr.mod)


# YUBR --------------------------------

# variable importance across the whole predictor set
yubr.varimp <- varimp.diag(y.data=as.numeric(yubr[,"flr"]), x.data=yubr[,xnames], ri.data=yubr[,"year"])

# variable selection
yubr.flr.mod.step <- bart.step(y.data=as.numeric(yubr[,"flr"]), x.data=yubr[,xnames], ri.data=yubr[,"year"], full=FALSE, quiet=TRUE) 
save(yubr.flr.mod.step, file="output/BART/bart.step.models.YUBR.Rdata")

# load(file="output/BART/bart.step.models.YUBR.Rdata")

summary(yubr.flr.mod.step)

varimp(yubr.flr.mod.step, plots=TRUE)
partial(yubr.flr.mod.step) # this needs debugging, per CJC
plot(yubr.flr.mod.step) #??

# refitting with vars indicated by varimp:
yubr.preds <- c("pptY0", "tminW0", "tmaxW0vW1", "tminW0vW1")

yubr.mod <- rbart_vi(as.formula(paste(paste('flr', paste(yubr.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yubr,
	group.by = yubr[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(yubr.mod)
varimp(yubr.mod)

save(yubr.mod, file="output/BART/bart.ri.model.YUBR.Rdata")


# YUJA --------------------------------

# variable importance across the whole predictor set
yuja.varimp <- varimp.diag(y.data=as.numeric(yuja[,"flr"]), x.data=yuja[,xnames], ri.data=yuja[,"year"])

# variable selection
yuja.flr.mod.step <- bart.step(y.data=as.numeric(yuja[,"flr"]), x.data=yuja[,xnames], ri.data=yuja[,"year"], full=FALSE, quiet=TRUE) 
save(yuja.flr.mod.step, file="output/BART/bart.step.models.YUJA.Rdata")

# load(file="output/BART/bart.step.models.YUJA.Rdata")

summary(yuja.flr.mod.step)

varimp(yuja.flr.mod.step)

# refitting with vars indicated by varimp:
yuja.preds <- c("vpdmaxW0vW1", "tmaxW0", "vpdmaxW0", "pptW0", "pptW0W1")


#------------------------------------------------
# AND ALSO


# plot the predictors in raw observations ...
jotr.preds <- c("tmaxW0vW1", "pptY0", "tmaxW0", "tminW0")

jotr.pred.plot <- flow2 %>% dplyr::select(year, flr, all_of(jotr.preds)) %>% pivot_longer(all_of(jotr.preds), names_to="Predictor", values_to="Value")


{cairo_pdf(file="output/figures/RImod_best_predictors_Jotr.pdf", width=6, height=2.5)

ggplot(jotr.pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.25, size=0.25) + geom_boxplot(alpha=0.5, aes(color=Predictor, fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,2,2,7)], guide="none") +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,2,2,7)], guide="none") +

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + dark_mode(theme_minimal()) + theme(plot.background=element_rect(color="black"))


}
dev.off()




# plot the predictors in raw observations ...
yubr.preds <- c("tmaxW0vW1", "pptY0", "vpdmaxW0vW1")

yubr.pred.plot <- yubr %>% dplyr::select(year, flr, all_of(yubr.preds)) %>% pivot_longer(all_of(yubr.preds), names_to="Predictor", values_to="Value")

{cairo_pdf(file="output/figures/RImod_best_predictors_YUBR.pdf", width=6, height=2.5)

ggplot(yubr.pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.25, size=0.25) + geom_boxplot(alpha=0.5, aes(color=Predictor, fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,2,7)], guide=FALSE) +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,2,7)], guide=FALSE) +

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + dark_mode(theme_minimal()) + theme(plot.background=element_rect(color="black"))

}
dev.off()




# plot the predictors in raw observations ...
yuja.preds <- c("pptW0", "vpdmaxW0", "vpdmaxW0vW1")

yuja.pred.plot <- yuja %>% dplyr::select(year, flr, all_of(yuja.preds)) %>% pivot_longer(all_of(yuja.preds), names_to="Predictor", values_to="Value")


{cairo_pdf(file="output/figures/RImod_best_predictors_YUJA.pdf", width=6, height=2.5)

ggplot(yuja.pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.25, size=0.25) + geom_boxplot(alpha=0.5, aes(color=Predictor, fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,7,7)], guide=FALSE) +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,7,7)], guide=FALSE) +

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + dark_mode(theme_minimal()) + theme(plot.background=element_rect(color="black"))


}
dev.off()


