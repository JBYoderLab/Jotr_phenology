# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2022.08.08

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("/Volumes/GoogleDrive/Other computers/My MacBook Pro 2020/Documents/Academic/Active_projects/Jotr_phenology")
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

# variant datasets -- dealing with the second flowering in 2019
flow2 <- flow %>% filter(year!=2019.5) # drop the weird observations
flow3 <- flow
flow3$year[flow3$year==2019.5] <- 2019 # or merge 2019.5 into 2019?

glimpse(flow2)

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="Eastern") 
yubr <- filter(flow2, type=="Western")

glimpse(yuja) # 1,048 obs
glimpse(yubr) # 1,280 obs

#-------------------------------------------------------------------------
# fit candidate BART models, stepwise

# predictors
xnames <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1", "vpdmaxW0", "vpdminW0", "vpdmaxW0vW1", "vpdminW0vW1") # year + weather data, "curated"

# flrmod <- bart(y.train=as.numeric(flow[,"flr"]), x.train=flow[,xnames], keeptrees=TRUE)

# summary(flrmod)


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


# YUJA --------------------------------

# variable importance across the whole predictor set
yuja.varimp <- varimp.diag(y.data=as.numeric(yuja[,"flr"]), x.data=yuja[,xnames], ri.data=yuja[,"year"])

# variable selection
yuja.flr.mod.step <- bart.step(y.data=as.numeric(yuja[,"flr"]), x.data=yuja[,xnames], ri.data=yuja[,"year"], full=FALSE, quiet=TRUE) 
save(yuja.flr.mod.step, file="output/BART/bart.step.models.YUJA.Rdata")

# load(file="output/BART/bart.step.models.YUJA.Rdata")

summary(yuja.flr.mod.step)
varimp(yuja.flr.mod.step)



#------------------------------------------------
# AND ALSO

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


