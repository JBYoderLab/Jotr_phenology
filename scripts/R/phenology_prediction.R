# Using BARTs to predict historical Joshua tree flowering
# last used/modified jby, 2022.07.12

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("embarcadero")

MojExt <- c(-119, -112, 33, 38) # need this later!!


#-----------------------------------------------------------
# load the data
flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, gridded and annualized
flow2 <- flow %>% filter(year!=2019.5) # drop the weird observations

glimpse(flow2) # 2,328 observations, 9 candidate predictors

yubr <- filter(flow2, type=="Western")
glimpse(yubr) # 1,280 observations for Y brevifolia
yuja <- filter(flow2, type=="Eastern")
glimpse(yuja) # 1,048 for Y jaegeriana



#-----------------------------------------------------------
# build species predictions from historical PRISM (let's do 1900-2022)

# Jotr --------------------------------
if(!dir.exists("output/BART/predictions")) dir.create("output/BART/predictions")

# need to load in the model at some point?!
# load(file="output/BART/bart.step.models.jotr.Rdata")

fixed <- attr(jotr.flr.mod.step$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(jotr.flr.mod.step, preds[[attr(jotr.flr.mod.step$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}



# YUBR --------------------------------
dir.create("output/BART/predictions.YUBR")

load(file="output/BART/bart.step.models.YUBR.Rdata")

fixed <- attr(yubr.flr.mod.step$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yubr.flr.mod.step, preds[[attr(yubr.flr.mod.step$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions.YUBR/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}

# YUJA --------------------------------
dir.create("output/BART/predictions.YUJA")

load(file="output/BART/bart.step.models.YUJA.Rdata")

fixed <- attr(yuja.flr.mod.step$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yuja.flr.mod.step, preds[[attr(yuja.flr.mod.step$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions.YUJA/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}
