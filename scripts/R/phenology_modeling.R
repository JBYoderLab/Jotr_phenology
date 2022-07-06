# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2022.07.05

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("/Volumes/GoogleDrive/Other computers/My MacBook Pro 2020/Documents/Academic/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")
# setwd("~/Jotr_phenology-main")

library("tidyverse")

library("embarcadero")

#-----------------------------------------------------------
# initial file loading

# flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized

flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized and normed

dim(flow) # 40 columns, whee
glimpse(flow)

# variant datasets -- dealing with the second flowering in 2019
flow2 <- flow %>% filter(year!=2019.5) # drop the weird observations
flow3 <- flow
flow3$year[flow3$year==2019.5] <- 2019 # or merge 2019.5 into 2019?

#-------------------------------------------------------------------------
# fit a BART model

# predictors
xnames <- c("year", "ppty0q1", "tmaxy0q1", "tminy0q1", "ppty1q1", "tmaxy1q1", "tminy1q1", "ppty1q2", "tmaxy1q2", "tminy1q2", "ppty1q3", "tmaxy1q3", "tminy1q3", "ppty1q4", "tmaxy1q4", "tminy1q4", "ppty2q4", "tmaxy2q4", "tminy2q4") # year + quarterly climate data, "curated"

# flrmod <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,xnames], keeptrees=TRUE)

# summary(flrmod)

# variable selection
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
flr.mod.step <- bart.step(y.data=as.numeric(flow2[,"flr"]), x.data=flow2[,xnames], full=FALSE, quiet=TRUE) 
save(flr.mod.step, file="output/BART/bart.step.models.Rdata")

# load(file="output/BART/bart.step.models.Rdata")

summary(flr.mod.step)

varimp(flr.mod.step, plots=TRUE) # hrm ...
varimp(flr.mod.step, plots=FALSE)

partial(flr.mod.step, x.vars=attr(flr.mod.step$fit$data@x, "term.labels"))
# hrm ...


#-------------------------------------------------------------------------
# Inspect kept predictors' relationships with flowering

attr(flr.mod.step$fit$data@x, "term.labels")

keptvars <- attr(flr.mod.step$fit$data@x, "term.labels")[-1] # informed by modeling above

flow.ln <- flow %>% dplyr::select(c("flr", "year", all_of(keptvars))) %>% pivot_longer(all_of(keptvars), names_to="Predictor", values_to="Value")

ggplot(flow.ln, aes(x=flr, y=Value)) + geom_boxplot() + facet_grid(Predictor~year, scale="free_y")

#-------------------------------------------------------------------------
# Fit a random-effect model with kept predictors

FLR.ri <- rbart_vi(as.formula(paste(paste('flr', paste(attr(flr.mod.step$fit$data@x, "term.labels"),  collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = flow,
	group.by = flow[,'year'],
	n.chains = 1,
#	k = SDMstep$fit$model@node.prior@k,
	power = flr.mod.step$fit$model@tree.prior@power,
	base = flr.mod.step$fit$model@tree.prior@base,
	keepTrees = TRUE)

summary(FLR.ri)

save(FLR.ri, file="output/BART/bart.ri.model.Rdata")

load("output/BART/bart.ri.model.Rdata")

embarcadero:::plot.ri(FLR.ri) # that looks familiar!









