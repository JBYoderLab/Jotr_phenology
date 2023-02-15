# Using BARTs to predict historical Joshua tree flowering
# last used/modified jby, 2023.01.04

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/Jotr_phenology")

library("tidyverse")

library("embarcadero")

MojExt <- c(-119, -112, 33, 38) # need this later!!


#-----------------------------------------------------------
# load the data
flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, gridded and annualized
flow2 <- flow %>% filter(!(year==2019.5 & flr==TRUE)) %>% mutate(year=floor(year)) # drop the late-flowering anomaly

glimpse(flow2) # 2,600 observations, 9 candidate predictors

yubr <- filter(flow2, type=="YUBR")
glimpse(yubr) # 1,440 observations for Y brevifolia
yuja <- filter(flow2, type=="YUJA")
glimpse(yuja) # 1,160 for Y jaegeriana



#-----------------------------------------------------------
# build species predictions from historical PRISM (let's do 1900-2022)

# Jotr --------------------------------
if(!dir.exists("output/BART/predictions")) dir.create("output/BART/predictions")

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

# need to load in the model at some point?!
# load(file="output/BART/bart.ri.model.Jotr.Rdata")

fixed <- attr(jotr.mod$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(jotr.mod, preds[[attr(jotr.mod$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}

# YUBR --------------------------------
dir.create("output/BART/predictions.YUBR")

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

fixed <- attr(yubr.mod$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yubr.mod, preds[[attr(yubr.mod$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions.YUBR/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}

# YUJA --------------------------------
dir.create("output/BART/predictions.YUJA")

# refitting with vars indicated by varimp:
yuja.preds <- c("vpdmaxW0vW1", "tmaxW0", "vpdmaxW0", "pptW0", "pptW0W1")

yuja.mod <- rbart_vi(as.formula(paste(paste('flr', paste(yuja.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yuja,
	group.by = yuja[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(yuja.mod)

fixed <- attr(yuja.mod$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yuja.mod, preds[[attr(yuja.mod$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions.YUJA/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}
