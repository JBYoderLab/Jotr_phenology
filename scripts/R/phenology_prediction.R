# Using BARTs to predict historical Joshua tree flowering
# last used/modified jby, 2023.03.11

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/Jotr_phenology")

library("tidyverse")

library("embarcadero")

MojExt <- c(-119, -112, 33, 38) # need this later!!


#-----------------------------------------------------------
# load the data
flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, gridded and annualized
flow2 <- flow %>% filter(!(year==2019.5 & flr==TRUE), year>=2008) %>% mutate(year=floor(year)) # drop the late-flowering anomaly

glimpse(flow2) # 3,016 observations, 13 candidate predictors

yubr <- filter(flow2, type=="YUBR")
glimpse(yubr) # 1,761 observations for Y brevifolia
yuja <- filter(flow2, type=="YUJA")
glimpse(yuja) # 1,335 for Y jaegeriana



#-----------------------------------------------------------
# build species predictions from historical PRISM (let's do 1900-2022)

# Jotr, no RI --------------------------------
if(!dir.exists("output/BART/predictions")) dir.create("output/BART/predictions")

# load the saved model developed in `phenology_modeling.R`
jotr.mod <- read_rds(file="output/BART/bart.model.Jotr.rds") # this now WORKS

jotr.preds <- attr(jotr.mod$fit$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(jotr.mod, preds[[jotr.preds]], splitby=20)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions/BART_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

} # END loop over years


# Jotr, with RI --------------------------------
if(!dir.exists("output/BART/RI.predictions")) dir.create("output/BART/RI.predictions")

# load the saved model developed in `phenology_modeling.R`
jotr.RImod <- read_rds(file="output/BART/bart.ri.model.Jotr.rds") # this now WORKS

jotr.preds <- attr(jotr.RImod$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(jotr.RImod, preds[[attr(jotr.RImod$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/RI.predictions/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

} # END loop over years



