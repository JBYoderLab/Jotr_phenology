# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2022.07.12

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("/Volumes/GoogleDrive/Other computers/My MacBook Pro 2020/Documents/Academic/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")
# setwd("~/Jotr_phenology-main")

library("tidyverse")

library("embarcadero")

#-----------------------------------------------------------
# initial file loading

# flow <- read.csv("output/flowering_obs_climate.csv") # flowering/not flowering, gridded and annualized
# flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized and normed

flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)

# variant datasets -- dealing with the second flowering in 2019
flow2 <- flow %>% filter(year!=2019.5) # drop the weird observations
flow3 <- flow
flow3$year[flow3$year==2019.5] <- 2019 # or merge 2019.5 into 2019?

glimpse(flow2)

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="Eastern") 
yubr <- filter(flow2, type=="Western")

glimpse(yuja)
glimpse(yubr)

#-------------------------------------------------------------------------
# fit candidate BART models, stepwise

# predictors
xnames <- c("year", "pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1") # year + weather data, "curated"

# flrmod <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,xnames], keeptrees=TRUE)

# summary(flrmod)

# variable selection
# YUBR --------------------------------
yubr.flr.mod.step <- bart.step(y.data=as.numeric(yubr[,"flr"]), x.data=yubr[,xnames], full=FALSE, quiet=TRUE) 
save(yubr.flr.mod.step, file="output/BART/bart.step.models.YUBR.Rdata")

# load(file="output/BART/bart.step.models.Rdata")

summary(yubr.flr.mod.step)
varimp(yubr.flr.mod.step, plots=FALSE)

partial(flr.mod.step, x.vars=attr(flr.mod.step$fit$data@x, "term.labels"))
# hrm ...

# YUJA --------------------------------
yuja.flr.mod.step <- bart.step(y.data=as.numeric(yuja[,"flr"]), x.data=yuja[,xnames], full=FALSE, quiet=TRUE) 
save(yuja.flr.mod.step, file="output/BART/bart.step.models.YUJA.Rdata")

# load(file="output/BART/bart.step.models.Rdata")

summary(yuja.flr.mod.step)
varimp(yuja.flr.mod.step, plots=FALSE)


#-------------------------------------------------------------------------
# Fit random-intercept models with kept predictors

# YUBR --------------------------------
yubr.FLR.ri <- rbart_vi(as.formula(paste(paste('flr', paste(attr(yubr.flr.mod.step$fit$data@x, "term.labels"),  collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yubr,
	group.by = yubr[,'year'],
	n.chains = 1,
#	k = SDMstep$fit$model@node.prior@k,
	power = yubr.flr.mod.step$fit$model@tree.prior@power,
	base = yubr.flr.mod.step$fit$model@tree.prior@base,
	keepTrees = TRUE)

summary(yubr.FLR.ri)

save(yubr.FLR.ri, file="output/BART/bart.ri.model.YUBR.Rdata")

load("output/BART/bart.ri.model.YUBR.Rdata")

embarcadero:::plot.ri(yubr.FLR.ri)

# YUJA --------------------------------
yuja.FLR.ri <- rbart_vi(as.formula(paste(paste('flr', paste(attr(yuja.flr.mod.step$fit$data@x, "term.labels"),  collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yuja,
	group.by = yuja[,'year'],
	n.chains = 1,
#	k = SDMstep$fit$model@node.prior@k,
	power = yuja.flr.mod.step$fit$model@tree.prior@power,
	base = yuja.flr.mod.step$fit$model@tree.prior@base,
	keepTrees = TRUE)

summary(yuja.FLR.ri)

save(yuja.FLR.ri, file="output/BART/bart.ri.model.YUJA.Rdata")

load("output/BART/bart.ri.model.YUJA.Rdata")

embarcadero:::plot.ri(yuja.FLR.ri)















