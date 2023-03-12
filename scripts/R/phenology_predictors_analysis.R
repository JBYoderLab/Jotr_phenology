# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2023.02.24

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
load("output/BART/bart.ri.varimp.Jotr.Rdata")


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

#-------------------------------------------------------------------------
# plot raw predictors and time-series and maybe also partials?

# plot the predictors in raw observations ...
jotr.preds <- c("tmaxW0vW1", "vpdmaxW0vW1", "pptY0", "vpdminW0vW1", "tmaxW0", "pptW0") # should match above

jotr.pred.plot <- flow2 %>% dplyr::select(year, flr, all_of(jotr.preds)) |> pivot_longer(all_of(jotr.preds), names_to="Predictor", values_to="Value") |> mutate(Predictor=factor(Predictor,jotr.preds))

{cairo_pdf(file="output/figures/BART_best_predictors_Jotr.pdf", width=6, height=2.5)

ggplot(jotr.pred.plot, aes(x=flr, y=Value)) + geom_jitter(alpha=0.25, size=0.25) + geom_boxplot(alpha=0.5, aes(color=Predictor, fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(2,3,1,3,2,1)], guide="none") +
scale_fill_manual(values=park_palette("JoshuaTree")[c(2,3,1,3,2,1)], guide="none") +

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + theme_minimal(base_size=9) + theme(plot.background=element_rect(color="black"))


}
dev.off()

#-------------------------------------------------------------------------
# time-series of selected predictors

jotr.hist.flowering <- read.csv("output/historic_flowering_reconst_jotr.csv")

sample_sites <- jotr.hist.flowering |> filter(year==2022) |> sample_n(300) |> dplyr::select(lat, lon)

sample_predictor_history <- data.frame(matrix(0,0,3+length(jotr.preds)))
names(sample_predictor_history) <- c("lat", "lon", "year", jotr.preds)

for(yr in 1900:2022){ # loop over years

sample_predictor_history <- rbind(sample_predictor_history, data.frame(sample_sites,year=yr, raster::extract(brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri",sep="")), sample_sites[,c("lon","lat")])[,jotr.preds]))

}

# now compile trends ...
pred_history_long <- sample_predictor_history |> pivot_longer(all_of(jotr.preds), names_to="predictor", values_to="value")

pred_history_long$predictor[pred_history_long$predictor=="pptY0"] <- "Year-of precip, mm"
pred_history_long$predictor[pred_history_long$predictor=="pptW0"] <- "Winter-of precip, mm"
pred_history_long$predictor[pred_history_long$predictor=="tmaxW0"] <- "Winter-of max temp, °C"
pred_history_long$predictor[pred_history_long$predictor=="tmaxW0vW1"] <- "Max temp winter-of vs winter before, °C"
pred_history_long$predictor[pred_history_long$predictor=="vpdmaxW0vW1"] <- "Max VPD winter-of vs winter before"
pred_history_long$predictor[pred_history_long$predictor=="vpdminW0vW1"] <- "Min VPD winter-of vs winter before"

pred_history_summary <- pred_history_long |> group_by(year, predictor) |> summarize(mn=mean(value), md=median(value), lo95=quantile(value, 0.025), up95=quantile(value, 0.975), lo50=quantile(value, 0.25), up50=quantile(value, 0.75))


{cairo_pdf(file="output/figures/BART_predictors_sampled_history.pdf", width=6, height=4)

ggplot() + geom_linerange(data=pred_history_summary, aes(x=year, ymin=lo95, ymax=up95), size=0.25, color="gray40") + geom_linerange(data=pred_history_summary, aes(x=year, ymin=lo50, ymax=up50), size=1.25, color="gray40") + geom_smooth(data=pred_history_long, aes(x=year, y=value), method="lm", size=0.5, color="blue") + geom_point(data=pred_history_summary, aes(x=year, y=md), color="white", size=0.25) + facet_wrap("predictor", scale="free_y") + labs(x="Year", y=NULL, title="Trends in best predictors, 1900-2022", subtitle="300 random points, 95% density, 50% density, and median") + theme_minimal(base_size=9)

}
dev.off()

{cairo_pdf(file="output/figures/RImod_predictors_sampled_history_trends.pdf", width=6, height=4)

ggplot(pred_history_long, aes(x=year, y=value)) + geom_smooth(method="lm", size=0.25, color="white") + facet_wrap("predictor", scale="free_y") + labs(x="Year", y=NULL, title="Trends in best predictors, 1900-2022", subtitle="300 random points") + theme_minimal(base_size=9)

}
dev.off()


