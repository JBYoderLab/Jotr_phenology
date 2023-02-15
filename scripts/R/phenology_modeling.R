# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2023.01.04

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/Jotr_phenology")
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
flow2 <- flow %>% filter(!(year==2019.5 & flr==TRUE)) %>% mutate(year=floor(year)) # drop the late-flowering anomaly
flow3 <- flow
flow3$year[flow3$year==2019.5] <- 2019 # or merge 2019.5 into 2019?

glimpse(flow2) # 2,600 in our final working set

ggplot(flow2, aes(x=lon, y=lat, color=flr)) + geom_point() + facet_wrap("year") + theme_bw()

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="YUJA") 
yubr <- filter(flow2, type=="YUBR")

glimpse(yuja) # 1,160 obs (after iffy ones excluded)
glimpse(yubr) # 1,440 obs

#-------------------------------------------------------------------------
# predictor exploration

xnames <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1", "vpdmaxW0", "vpdminW0", "vpdmaxW0vW1", "vpdminW0vW1") # year + weather data, "curated"

flow2.ln <- flow2 %>% pivot_longer(cols=all_of(xnames), values_to="value", names_to="variable")
glimpse(flow2.ln)

{cairo_pdf("output/figures/predictor_differences.pdf", width=11, height=4)

ggplot(flow2.ln, aes(x=value, y=1+flr)) + geom_jitter(alpha=0.05, height=0.125) + geom_boxplot(alpha=0.1, aes(group=flr)) + geom_smooth(method="gam") + facet_wrap("variable", scale="free_x", nrow=2) + scale_y_continuous(breaks=1:2, labels=c("False", "True")) + labs(x="Variable value", y="Flowers observed?") + theme_bw()

}
dev.off()

#-------------------------------------------------------------------------
# fit candidate BART models, stepwise

# predictors
xnames <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1", "vpdmaxW0", "vpdminW0", "vpdmaxW0vW1", "vpdminW0vW1") # year + weather data, "curated"

# Full range ------------------------------------

# variable importance across the whole predictor set
jotr.varimp <- varimp.diag(y.data=as.numeric(flow2[,"flr"]), x.data=flow2[,xnames], ri.data=flow2[,"year"])

save(jotr.varimp, file="output/BART/bart.ri.varimp.Jotr.Rdata")


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

save(jotr.mod, file="output/BART/bart.ri.model.Jotr.Rdata")

# load(file="output/BART/bart.ri.model.Jotr.Rdata")

summary(jotr.mod)


# YUBR --------------------------------

# fitting with vars indicated by varimp on global data:
yubr.preds <- jotr.preds

yubr.mod <- rbart_vi(as.formula(paste(paste('flr', paste(yubr.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yubr,
	group.by = yubr[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(yubr.mod)
varimp(yubr.mod)

save(yubr.mod, file="output/BART/bart.ri.model.YUBR.Rdata")


# YUJA --------------------------------

# refitting with vars indicated by varimp:
yuja.preds <- jotr.preds

yuja.mod <- rbart_vi(as.formula(paste(paste('flr', paste(yuja.preds, collapse=' + '), sep = ' ~ '), 'year', sep=' - ')),
	data = yuja,
	group.by = yuja[,'year'],
	n.chains = 1,
	k = 2,
	power = 2,
	base = 0.95,
	keepTrees = TRUE)

summary(yuja.mod)
varimp(yuja.mod)

save(yuja.mod, file="output/BART/bart.ri.model.YUJA.Rdata")

#------------------------------------------------
# AND ALSO


# plot the predictors in raw observations ...
jotr.preds <- c("tmaxW0vW1", "pptY0", "tmaxW0", "tminW0")

jotr.pred.plot <- flow2 %>% dplyr::select(year, flr, all_of(jotr.preds)) %>% pivot_longer(all_of(jotr.preds), names_to="Predictor", values_to="Value")


{cairo_pdf(file="output/figures/RImod_best_predictors_Jotr.pdf", width=6, height=2.5)

ggplot(jotr.pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.25, size=0.25) + geom_boxplot(alpha=0.5, aes(color=Predictor, fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,2,2,7)], guide="none") +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,2,2,7)], guide="none") +

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + dark_mode(theme_minimal()) + theme(plot.background=element_rect(color="black"))


}
dev.off()




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


