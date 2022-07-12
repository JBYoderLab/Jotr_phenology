# Using BARTs to predict historical Joshua tree flowering
# last used/modified jby, 2022.07.12

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("embarcadero")

MojExt <- c(-119, -112, 33, 38) # need this later!!


#-----------------------------------------------------------
# fit the RI model

# load the data
flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, gridded and annualized
flow2 <- flow %>% filter(year!=2019.5) # drop the weird observations

glimpse(flow2) # 2,328 observations, 9 candidate predictors

yubr <- filter(flow2, type=="Western")
glimpse(yubr) # 1,280 observations for Y brevifolia
yuja <- filter(flow2, type=="Eastern")
glimpse(yuja) # 1,048 for Y jaegeriana


# load the step-fitted models
load(file="output/BART/bart.step.models.YUBR.Rdata")
load(file="output/BART/bart.step.models.YUJA.Rdata")

# YUBR --------------------------------
yubr.preds <- attr(yubr.flr.mod.step$fit$data@x, "term.labels")
if(!"year" %in% yubr.preds) yubr.preds <- c("year", yubr.preds)

# fit the RI model with year as RI
yubr.FLR.ri <- rbart_vi(as.formula(paste(paste('flr', paste(yubr.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yubr,
	group.by = yubr[,'year'],
	n.chains = 1,
#	k = SDMstep$fit$model@node.prior@k,
	power = yubr.flr.mod.step$fit$model@tree.prior@power,
	base = yubr.flr.mod.step$fit$model@tree.prior@base,
	keepTrees = TRUE)

summary(yubr.FLR.ri)

save(yubr.FLR.ri, file="output/BART/bart.ri.model.YUBR.Rdata")

# YUJA --------------------------------
yuja.preds <- attr(yuja.flr.mod.step$fit$data@x, "term.labels")
if(!"year" %in% yuja.preds) yuja.preds <- c("year", yuja.preds)

# fit the RI model with year as RI
yuja.FLR.ri <- rbart_vi(as.formula(paste(paste('flr', paste(yuja.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yuja,
	group.by = yuja[,'year'],
	n.chains = 1,
#	k = SDMstep$fit$model@node.prior@k,
	power = yuja.flr.mod.step$fit$model@tree.prior@power,
	base = yuja.flr.mod.step$fit$model@tree.prior@base,
	keepTrees = TRUE)

summary(yuja.FLR.ri)

save(yuja.FLR.ri, file="output/BART/bart.ri.model.YUJA.Rdata")


#-----------------------------------------------------------
# build species predictions from historical PRISM (let's do 1900-2022)

# YUBR --------------------------------
dir.create("output/BART/predictions.YUBR")

fixed <- attr(yubr.FLR.ri$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

pptW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep="")), raster)),sum)
pptY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:2,".bil", sep="")), raster)),sum)
pptW0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0Y1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q",4:2,".bil", sep="")), raster)),sum)
tmaxW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q4.bil", sep="")), raster)), max)
tminW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q4.bil", sep="")), raster)), min)
tmaxW0vW1 <- tmaxW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-2,"Q4.bil", sep="")), raster)), max)
tminW0vW1 <- tminW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-2,"Q4.bil", sep="")), raster)), min)

preds <- brick(c(pptW0, pptY0, pptW0W1, pptY0W1, pptY0Y1, tmaxW0, tminW0, tmaxW0vW1, tminW0vW1))
names(preds) <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1")

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yubr.FLR.ri, preds[[attr(yubr.FLR.ri$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions.YUBR/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}

# YUJA --------------------------------
dir.create("output/BART/predictions.YUJA")

fixed <- attr(yuja.FLR.ri$fit[[1]]$data@x, "term.labels")

# LOOP over years
for(yr in 1900:2022){

# yr <- 1900

pptW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep="")), raster)),sum)
pptY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:2,".bil", sep="")), raster)),sum)
pptW0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0Y1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q",4:2,".bil", sep="")), raster)),sum)
tmaxW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q4.bil", sep="")), raster)), max)
tminW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q4.bil", sep="")), raster)), min)
tmaxW0vW1 <- tmaxW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-2,"Q4.bil", sep="")), raster)), max)
tminW0vW1 <- tminW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-2,"Q4.bil", sep="")), raster)), min)

preds <- brick(c(pptW0, pptY0, pptW0W1, pptY0W1, pptY0Y1, tmaxW0, tminW0, tmaxW0vW1, tminW0vW1))
names(preds) <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1")

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(yuja.FLR.ri, preds[[attr(yuja.FLR.ri$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions.YUJA/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}
