# Using BARTs to predict historical Joshua tree flowering
# last used/modified jby, 2022.07.05

# rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("embarcadero")

#-----------------------------------------------------------
# fit the RI model

MojExt <- c(-119, -112, 33, 38) # need this later!!

# load the data
flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized
flow2 <- flow %>% filter(year!=2019.5) # drop the weird observations

dim(flow2) # 40 columns, 2,063 observations, whee
glimpse(flow2)


# load the step-fitted model
load(file="output/BART/bart.step.models.Rdata")

summary(flr.mod.step) # huh, no year

# partial(flr.mod.step.v2, x.vars=c("year", "ppty1q1", "ppty1q3", "ppty2q1", "tmaxy2q4"))

preds <- attr(flr.mod.step$fit$data@x, "term.labels")
if(!"year" %in% preds) preds <- c("year", preds)

# fit the RI model with year as RI
FLR.ri <- rbart_vi(as.formula(paste(paste('flr', paste(preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = flow,
	group.by = flow[,'year'],
	n.chains = 1,
#	k = SDMstep$fit$model@node.prior@k,
	power = flr.mod.step$fit$model@tree.prior@power,
	base = flr.mod.step$fit$model@tree.prior@base,
	keepTrees = TRUE)

summary(FLR.ri)

save(FLR.ri, file="output/BART/bart.ri.model.Rdata")

embarcadero:::plot.ri(FLR.ri) # that looks familiar!


#-----------------------------------------------------------
# LOOP over the full period (let's do 1901-2021)

fixed <- attr(FLR.ri$fit[[1]]$data@x, "term.labels")

dir.create("output/BART/predictions")

# LOOP over years
for(yr in 1901:2021){

# yr <- 2022

# parse predictors into useful variables
pds <- data.frame(v=gsub("(\\w+)y.+", "\\1", fixed), o=as.numeric(gsub(".+y(\\d).+", "\\1", fixed)), q=gsub(".+q(\\d)", "\\1", fixed))

# get list of PRISM data files for those variables
files <- paste("data/PRISM/annual_normed_1981-2010/", pds$v, "_normed_Mojave_", yr-pds$o, "Q", pds$q, ".bil", sep="")

# assemble a RasterStack with those data
predstack <- raster::stack(sapply(files, function(x) crop(raster::raster(x), MojExt)))
names(predstack) <- fixed

# prediction with the RI predictor (year) removed
pred.ri0 <- predict(FLR.ri, predstack[[attr(FLR.ri$fit[[1]]$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

pred.ri0 # BOOM

writeRaster(pred.ri0, paste("output/BART/predictions/BART_RI_predicted_flowering_",yr,".bil", sep=""), overwrite=TRUE)

}

