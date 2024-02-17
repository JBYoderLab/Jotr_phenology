# Using BARTs to predict historical Joshua tree flowering
# last used/modified jby, 2024.01.22

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/Jotr_phenology")

library("tidyverse")

library("embarcadero")

MojExt <- c(-119, -112, 33, 40) # need this later!!


#-----------------------------------------------------------
# load the data
flow <- read.csv("output/flowering_obs_climate_subsp.csv") # flowering/not flowering, gridded and annualized
flow2 <- flow %>% filter(!(year==2018.5 & flr==TRUE), year>=2008) %>% mutate(year=floor(year)) # drop the late-flowering anomaly

glimpse(flow2) # 3,148 records with 2023 now
filter(flow2, year!=2023) %>% glimpse() # 2,632 observations, 13 candidate predictors

yubr <- filter(flow2, type=="YUBR")
glimpse(yubr) 
filter(yubr, year!=2023) %>% glimpse() # 1,460 observations for Y brevifolia

yuja <- filter(flow2, type=="YUJA")
glimpse(yuja) 
filter(yuja, year!=2023) %>% glimpse() # 1,172 for Y jaegeriana


#-----------------------------------------------------------
# build species predictions from historical PRISM (let's do 1900-2022)

# Jotr, no RI --------------------------------
if(!dir.exists("output/BART/predictions")) dir.create("output/BART/predictions")

# load the saved model developed in `phenology_modeling.R`
jotr.mod <- read_rds(file="output/BART/bart.model.Jotr.rds") # this now WORKS

summary(jotr.mod)

jotr.preds <- attr(jotr.mod$fit$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2023){

# yr <- 2023

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
for(yr in 1900:2023){

# yr <- 2023

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(jotr.RImod, preds[[attr(jotr.RImod$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/RI.predictions/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

} # END loop over years


# YUJA, no RI --------------------------------
if(!dir.exists("output/BART/yuja.predictions")) dir.create("output/BART/yuja.predictions")

# load the saved model developed in `phenology_modeling.R`
yuja.mod <- read_rds(file="output/BART/bart.model.yuja.rds") # this now WORKS

summary(yuja.mod)

yuja.preds <- attr(yuja.mod$fit$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2023){

# yr <- 2023

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yuja.mod, preds[[yuja.preds]], splitby=20)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/yuja.predictions/BART_yuja_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

} # END loop over years

# YUBR, no RI --------------------------------
if(!dir.exists("output/BART/yubr.predictions")) dir.create("output/BART/yubr.predictions")

# load the saved model developed in `phenology_modeling.R`
yubr.mod <- read_rds(file="output/BART/bart.model.yubr.rds") # this now WORKS

summary(yubr.mod)

yubr.preds <- attr(yubr.mod$fit$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2023){

# yr <- 2023

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yubr.mod, preds[[yubr.preds]], splitby=20)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/yubr.predictions/BART_yubr_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

} # END loop over years
