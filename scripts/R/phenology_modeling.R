# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2023.03.31

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Jotr_phenology-main")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("embarcadero")
library("ggdark")
library("raster")

source("../shared/Rscripts/base_graphics.R")

#-----------------------------------------------------------
# initial file loading

# flow <- read.csv("output/flowering_obs_climate.csv") # flowering/not flowering, gridded and annualized
# flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized and normed

flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)


# variant datasets -- dealing with the second flowering in 2018
flow2 <- flow %>% filter(!(year==2018.5 & flr==TRUE), year>=2008) %>% mutate(year=floor(year)) # drop the late-flowering anomaly
flow3 <- flow |> filter(year>=2008)
flow3$year[flow3$year==2018.5] <- 2018 # or merge 2018.5 into 2018?

glimpse(flow2) # 3,009 in our final working set

table(flow2$year)

flow4 <- flow2 |> filter(year>=2016) # years with at least 100 records
glimpse(flow4) # 2,674 from 2016 on

# ggplot(flow2, aes(x=lon, y=lat, color=flr)) + geom_point() + facet_wrap("year") + theme_bw()

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="YUJA") 
yubr <- filter(flow2, type=="YUBR")

glimpse(yuja) # 1,310 obs (after iffy ones excluded)
glimpse(yubr) # 1,699 obs

#-------------------------------------------------------------------------
# predictor exploration

xnames <- c("pptY0", "pptY1", "pptY2", "pptY0Y1", "pptY1Y2", "tmaxY0", "tmaxY0Y1", "tminY0", "tminY0Y1", "vpdmaxY0", "vpdmaxY0Y1", "vpdminY0", "vpdminY0Y1") # weather data, curated

flow2.ln <- flow2 %>% pivot_longer(cols=all_of(xnames), values_to="value", names_to="variable") |> mutate(variable = factor(variable, xnames))
glimpse(flow2.ln)

{cairo_pdf("output/figures/predictor_differences.pdf", width=11, height=4)

ggplot(flow2.ln, aes(x=value, y=1+flr)) + geom_jitter(alpha=0.05, height=0.125) + geom_boxplot(alpha=0.1, aes(group=flr)) + geom_smooth(method="gam") + facet_wrap("variable", scale="free_x", nrow=2) + scale_y_continuous(breaks=1:2, labels=c("False", "True")) + labs(x="Variable value", y="Flowers observed?") + theme_bw()

}
dev.off()

#-------------------------------------------------------------------------
# fit candidate BART models, stepwise

# predictors
xnames <- c("pptY0", "pptY1", "pptY2", "pptY0Y1", "pptY1Y2", "tmaxY0", "tmaxY0Y1", "tminY0", "tminY0Y1", "vpdmaxY0", "vpdmaxY0Y1", "vpdminY0", "vpdminY0Y1") # weather data, curated

# Full range ------------------------------------

# variable importance across the whole predictor set
jotr.varimp <- varimp.diag(y.data=as.numeric(flow2[,"flr"]), x.data=flow2[,xnames]) # now sans RI
# favors pptY1Y2, pptY0Y1, vpdmaxY0, vpdminY0Y1, tminY0, tmaxY0Y1

write_rds(jotr.varimp, file="output/BART/bart.varimp.Jotr.rds") # switching save modes now
# jotr.varimp <- read_rds("output/BART/bart.varimp.Jotr.rds")

# and, sure, let's do stepwise fitting too
jotr.mod.step <- bart.step(y.data=as.numeric(flow2[,"flr"]), x.data=flow2[,xnames], full=FALSE, quiet=TRUE)
# favors pptY0Y1 pptY1Y2 tmaxY0Y1 tminY0 tminY0Y1 vpdmaxY0 vpdminY0Y1 (holy shit, that's all the varimp calls?)

invisible(jotr.mod.step$fit$state)
write_rds(jotr.mod.step, file="output/BART/bart.step.model.Jotr.rds")

summary(jotr.mod.step)

# refitting with vars indicated by varimp:
jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

jotr.mod <- bart(y.train=as.numeric(flow2[,"flr"]), x.train=flow2[,jotr.preds], keeptrees=TRUE)

summary(jotr.mod)

invisible(jotr.mod$fit$state)
write_rds(jotr.mod, file="output/BART/bart.model.Jotr.rds")
# jotr.mod <- read_rds("output/BART/bart.model.Jotr.rds")

p <- partial(jotr.mod, jotr.preds, trace=FALSE) # visualize partials
varimp(jotr.mod)

pd2_36 <- pd2bart(y.train=as.numeric(flow2[,"flr"]), x.train=flow2[,jotr.preds], xind=c(3,6))  # frowning


# now, spartials ...

for(yr in unique(flow2$year)){

# yr <- 2022

{cairo_pdf(paste("output/figures/BART_spartials_jotr_", yr, ".pdf", sep=""), width=9, height=6.5)
plot(spartial(jotr.mod, brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_", yr, ".grd", sep="")), x.vars=jotr.preds))
}
dev.off()

}


#-------------------------------------------------------------------------
# LOO by year, for more confirmation

table(flow2$year) # do we have enough for all years? Yeah sure

LOOvalid <- data.frame(matrix(0,0,4))
colnames(LOOvalid) <- c("year", "N_flr", "N_noflr", "AUC")

# LOOP over years
for(yr in unique(flow2$year)){

# yr <- 2010

inbag <- flow2 |> filter(year!=yr)
oobag <- flow2 |> filter(year==yr)

testmod <- bart(y.train=as.numeric(inbag[,"flr"]), x.train=inbag[,jotr.preds], keeptrees=TRUE)

# data for OOB year
OOBpreds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction at OOB sites with the RI predictor (year) removed
testpred.ri0 <- predict(testmod, OOBpreds[[attr(testmod$fit$data@x, "term.labels")]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

OOBpreds <- raster::extract(testpred.ri0, oobag[,c("lon", "lat")]) # predicted OOB sites with model

# hacked out of the summary function for rbarts
auc <- performance(prediction(OOBpreds, oobag$flr),"auc")@y.values[[1]]

LOOvalid <- rbind(LOOvalid, data.frame(year=yr, N_flr=length(which(oobag$flr)), N_noflr=length(which(!oobag$flr)), AUC=auc))

} # END loop over years

LOOvalid <- LOOvalid |> arrange(year) # eeeeeh

write.table(LOOvalid, "output/BART/year-year-LOO.csv", sep=",", col.names=TRUE, row.names=FALSE)

mean(LOOvalid$AUC) # 0.59
sd(LOOvalid$AUC) # 0.17
sd(LOOvalid$AUC)/sqrt(nrow(LOOvalid)) # SE = 0.04

#-------------------------------------------------------------------------
# "recent" model based only on years with >100 observations

flow4 <- flow2 |> filter(year >= 2016)

jotr.recent.mod <- bart(y.train=as.numeric(flow4[,"flr"]), x.train=flow4[,jotr.preds], keeptrees=TRUE)

summary(jotr.recent.mod)

invisible(jotr.recent.mod$fit$state)
write_rds(jotr.recent.mod, file="output/BART/bart.recent.model.Jotr.rds")
# jotr.mod <- read_rds("output/BART/bart.recent.model.Jotr.rds")

# cross-validate with earlier data ---------------
earlyValid <- data.frame(matrix(0,0,4))
colnames(earlyValid) <- c("year", "lat", "lon", "flr", "PrFlr")

# LOOP over years
for(yr in 2008:2015){

# yr <- 2010
oobag <- flow2 |> filter(year==yr)

# data for OOB year
OOBpreds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction at OOB sites with the RI predictor (year) removed
testpred <- predict(jotr.recent.mod, OOBpreds[[attr(jotr.recent.mod$fit$data@x, "term.labels")]])

earlyPreds <- raster::extract(testpred, oobag[,c("lon", "lat")]) # predicted OOB sites with model

earlyValid <- rbind(earlyValid, data.frame(year=yr, lon=oobag$lon, lat=oobag$lat, flr=oobag$flr, PrFlr=earlyPreds[,"layer"]))

}

glimpse(earlyValid)

# hacked out of the summary function for rbarts
auc <- performance(prediction(earlyValid$PrFlr, earlyValid$flr), "auc")@y.values[[1]]

auc # 0.63, okay?


#-------------------------------------------------------------------------
# finally, now fit RI model
jotr.RImod <- rbart_vi(as.formula(paste(paste('flr', paste(jotr.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = flow2,
	group.by = flow2[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(jotr.RImod)

invisible(jotr.RImod$fit[[1]]$state) # MUST do this to save
write_rds(jotr.RImod, file="output/BART/bart.ri.model.Jotr.rds") # write out for downstream use


# jotr.RImod <-  read_rds(file="output/BART/bart.ri.model.Jotr.rds")

summary(jotr.RImod)


# LOO validation by year
LOOvalidRI <- data.frame(matrix(0,0,4))
colnames(LOOvalidRI) <- c("year", "N_flr", "N_noflr", "AUC")

# LOOP over years
for(yr in unique(flow2$year)){

# yr <- 2010

inbag <- flow2 |> filter(year!=yr)
oobag <- flow2 |> filter(year==yr)

testmod <- rbart_vi(as.formula(paste(paste('flr', paste(jotr.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = inbag,
	group.by = inbag[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

# data for OOB year
OOBpreds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri", sep=""))

# prediction at OOB sites with the RI predictor (year) removed
testpred.ri0 <- predict(testmod, OOBpreds[[jotr.preds]], splitby=20, ri.data=yr, ri.name='year', ri.pred=FALSE)

OOBpreds <- raster::extract(testpred.ri0, oobag[,c("lon", "lat")]) # predicted OOB sites with model

# hacked out of the summary function for rbarts
auc <- performance(prediction(OOBpreds, oobag$flr),"auc")@y.values[[1]]

LOOvalidRI <- rbind(LOOvalidRI, data.frame(year=yr, N_flr=length(which(oobag$flr)), N_noflr=length(which(!oobag$flr)), AUC=auc))

} # END loop over years

LOOvalidRI <- LOOvalidRI |> arrange(year) # eeeeeh

write.table(LOOvalidRI, "output/BART/RI_year-year-LOO.csv", sep=",", col.names=TRUE, row.names=FALSE)

mean(LOOvalidRI$AUC) # 0.57
sd(LOOvalidRI$AUC) # 0.15
sd(LOOvalidRI$AUC)/sqrt(nrow(LOOvalidRI)) # SE = 0.04




