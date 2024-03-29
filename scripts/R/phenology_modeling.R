# Using BARTs to model Joshua tree flowering
# best run on MAJEL
# last used/modified jby, 2023.12.22

rm(list=ls())  # Clears memory of all objects -- useful for debugging! But doesn't kill packages.

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Jotr_phenology-main")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("embarcadero")
library("ggdark")
library("raster")
library("sf")

source("../shared/Rscripts/base_graphics.R")

#-----------------------------------------------------------
# initial file loading

# flow <- read.csv("output/flowering_obs_climate.csv") # flowering/not flowering, gridded and annualized
# flow <- read.csv("output/flowering_obs_climate_normed.csv") # flowering/not flowering, gridded and annualized and normed

flow <- read.csv("output/flowering_obs_climate_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)


# variant datasets -- dealing with the second flowering in 2018
flow2 <- flow %>% filter(!(year==2018.5 & flr==TRUE), year>=2008) %>% mutate(year=floor(year)) # drop the late-flowering anomaly
flow3 <- flow |> filter(year>=2008)
flow3$year[flow3$year==2018.5] <- 2018 # or merge 2018.5 into 2018?

glimpse(flow2) # 2,632 in our final working set

table(flow2$year)

flow4 <- flow2 |> filter(year>=2016) # years with at least 100 records
glimpse(flow4) # 2,319 from 2016 on

# ggplot(flow2, aes(x=lon, y=lat, color=flr)) + geom_point() + facet_wrap("year") + theme_bw()

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="YUJA") 
yubr <- filter(flow2, type=="YUBR")

glimpse(yuja) # 1,172 obs (after iffy ones excluded)
glimpse(yubr) # 1,460 obs

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

# publication-ready figure assembly ...
glimpse(jotr.varimp$data)

jotr.varimp$data <- jotr.varimp$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

levels(jotr.varimp$data$names) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

jotr.varimp$labels$group <- "Trees"
jotr.varimp$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

{cairo_pdf(file="output/figures/varimp_Jotr.pdf", width=6, height=5)

jotr.varimp + scale_x_discrete(label=label_parse) + theme(legend.position=c(0.8, 0.7), axis.text=element_text(size=13), legend.text=element_text(size=12), legend.title=element_text(size=13))  # okay nice

}
dev.off()


# and, sure, let's do stepwise fitting too
jotr.mod.step <- bart.step(y.data=as.numeric(flow2[,"flr"]), x.data=flow2[,xnames], full=FALSE, quiet=TRUE)
# favors pptY1 pptY0Y1 pptY1Y2 tmaxY0Y1 tminY0 tminY0Y1 vpdmaxY0 vpdminY0Y1 (ooh, that's all the varimp calls?)

invisible(jotr.mod.step$fit$state)
write_rds(jotr.mod.step, file="output/BART/bart.step.model.Jotr.rds")
# jotr.mod.step <- read_rds("output/BART/bart.step.model.Jotr.rds")

summary(jotr.mod.step) # AUC is up, interesting

# refitting with vars indicated by varimp:
jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

jotr.mod <- bart(y.train=as.numeric(flow2[,"flr"]), x.train=flow2[,jotr.preds], keeptrees=TRUE)

summary(jotr.mod)

invisible(jotr.mod$fit$state)
write_rds(jotr.mod, file="output/BART/bart.model.Jotr.rds")
# jotr.mod <- read_rds("output/BART/bart.model.Jotr.rds")

p <- partial(jotr.mod, jotr.preds, trace=FALSE, smooth=5) # visualize partials
varimp(jotr.mod)

write_rds(p, file="output/BART/bart.model.Jotr.partials.rds")
# p <- read_rds("output/BART/bart.model.Jotr.partials.rds")

partplot <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT", p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT", p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]", p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD", p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]", p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp", p[[6]]$data)
				)


# now, spartials ...

for(yr in unique(flow2$year)){

# yr <- 2022

spYr <- spartial(jotr.mod, brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_", yr, ".grd", sep="")), x.vars=jotr.preds)

{cairo_pdf(paste("output/figures/BART_spartials_jotr_", yr, ".pdf", sep=""), width=9, height=6.5)
plot(spYr)
}
dev.off()

writeRaster(spYr, paste("output/BART/Jotr_BART_spartials_", yr,".grd", sep=""), overwrite=TRUE) # confirmed write-out and read-in
}

#-------------------------------------------------------------------------
# Partials and spartials in example years

# and then spartials ...
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range.shp")

states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_admin_1_states_provinces")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")


# 2019, "good year" example
goodSpart <- brick("output/BART/Jotr_BART_spartials_2019.grd")
projection(goodSpart)<-CRS("+init=epsg:4269")

goodSpart.mask <- mask(goodSpart, st_transform(sdm.pres[,2], crs=4269))

goodSpart.df <- cbind(coordinates(goodSpart.mask), as.data.frame(goodSpart.mask)) %>% filter(!is.na(pptY1Y2)) |> rename(lon=x, lat=y) |> pivot_longer(all_of(jotr.preds), names_to="predictor", values_to="prFL") |> mutate(predictor=factor(predictor,jotr.preds))

levels(goodSpart.df$predictor) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

goodex <- ggplot(goodSpart.df, aes(x=lon, y=lat, fill=prFL)) + geom_tile() + facet_wrap("predictor", nrow=3, labeller="label_parsed") + 
	scale_fill_gradient(low="#ffffcc", high="#253494", name="Marginal Pr(Flowers)", breaks=c(0.5,0.55,0.6,0.65), limits=c(0.5,0.675), labels=c("", 0.55, 0.6, 0.65)) + 
	labs(title="Spatial partial effects for 2019") + 
	theme_minimal() + theme(axis.title=element_blank(), axis.text=element_blank(), legend.position="bottom", legend.text=element_text(size=9), panel.background=element_rect(fill="white"), panel.grid=element_blank())

goodex

# map of observations
goodobs <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=states, fill="antiquewhite2", color="antiquewhite4") + 
	geom_sf(data=sdm.pres, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
	geom_point(data=filter(flow2, year==2019), aes(x=lon, y=lat, shape=flr, color=flr, size=flr), alpha=0.75) +
	scale_shape_manual(values=c(21,20), labels=c("No flowering", "Flowering"), name="Flowering") +
	scale_size_manual(values=c(1.1,1.3), labels=c("No flowering", "Flowering"), name="Flowering") +
	scale_color_manual(values=c("#ffffcc", "#253494"), labels=c("No flowering", "Flowering"), name="Flowering") + # nb these are false, true
	coord_sf(xlim = c(-119, -112.75), ylim = c(32.75, 38), expand = FALSE) +
	labs(title="2019 observations") +
	theme_bw(base_size=11) + theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "slategray3"), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position=c(0.5, 0.1), legend.box.margin=unit(c(0.005, 0.005, 0.005, 0.005), "in"), legend.title=element_blank(), legend.text=element_text(size=9), plot.margin=unit(c(0.05,0.05,0.05,0.05), "in"), legend.key.size=unit(0.075, "in"), legend.key=element_rect(fill="#A8A378"), legend.spacing=unit(c(0.005, 0.005, 0.005, 0.005), "in"), legend.direction="horizontal", legend.background=element_rect(color="black", linewidth=0.1))

# figure of good observations and spartials
{cairo_pdf("output/figures/present_goodex_spartials.pdf", width=10, height=5)

ggdraw() + draw_plot(goodobs, 0, 0.1, 0.4, 0.9) + draw_plot(goodex, 0.4, 0, 0.6, 0.99)

}
dev.off()


# 2020, "bad year" example
badSpart <- brick("output/BART/Jotr_BART_spartials_2020.grd")
projection(badSpart)<-CRS("+init=epsg:4269")

badSpart.mask <- mask(badSpart, st_transform(sdm.pres[,2], crs=4269))

badSpart.df <- cbind(coordinates(badSpart.mask), as.data.frame(badSpart.mask)) %>% filter(!is.na(pptY1Y2)) |> rename(lon=x, lat=y) |> pivot_longer(all_of(jotr.preds), names_to="predictor", values_to="prFL") |> mutate(predictor=factor(predictor,jotr.preds))

levels(badSpart.df$predictor) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

badex <- ggplot(badSpart.df, aes(x=lon, y=lat, fill=prFL)) + geom_tile() + facet_wrap("predictor", nrow=3, labeller="label_parsed") + scale_fill_gradient(low="#ffffcc", high="#253494", name="Marginal Pr(Flowers)", breaks=c(0.5,0.55,0.6,0.65), limits=c(0.5,0.675), labels=c("", 0.55, 0.6, 0.65)) + labs(title="Spatial partial effects for 2020") + theme_minimal() + theme(axis.title=element_blank(), axis.text=element_blank(), legend.position="bottom", legend.text=element_text(size=9), panel.background=element_rect(fill="white"), panel.grid=element_blank())

badex

# map of observations
badobs <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=states, fill="antiquewhite2", color="antiquewhite4") + 
	geom_sf(data=sdm.pres, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
	geom_point(data=filter(flow2, year==2020), aes(x=lon, y=lat, shape=flr, color=flr, size=flr), alpha=0.75) +
	scale_shape_manual(values=c(21,20), labels=c("No flowering", "Flowering"), name="Flowering") +
	scale_size_manual(values=c(1.1,1.3), labels=c("No flowering", "Flowering"), name="Flowering") +
	scale_color_manual(values=c("#ffffcc", "#253494"), labels=c("No flowering", "Flowering"), name="Flowering") + # nb these are false, true
	coord_sf(xlim = c(-119, -112.75), ylim = c(32.75, 38), expand = FALSE) +
	labs(title="2020 observations") +
	theme_bw(base_size=11) + theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "slategray3"), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position=c(0.5, 0.1), legend.box.margin=unit(c(0.005, 0.005, 0.005, 0.005), "in"), legend.title=element_blank(), legend.text=element_text(size=9), plot.margin=unit(c(0.05,0.05,0.05,0.05), "in"), legend.key.size=unit(0.075, "in"), legend.key=element_rect(fill="#A8A378"), legend.spacing=unit(c(0.005, 0.005, 0.005, 0.005), "in"), legend.direction="horizontal", legend.background=element_rect(color="black", linewidth=0.1))

badobs

# figure of bad-year observations and spartials
{cairo_pdf("output/figures/present_badex_spartials.pdf", width=10, height=5)

ggdraw() + draw_plot(badobs, 0, 0.1, 0.4, 0.9) + draw_plot(badex, 0.4, 0, 0.6, 0.99)

}
dev.off()




# lower panel assembly
spartials <- goodex + badex + plot_layout(guides="collect") & theme(legend.position="bottom", legend.text=element_text(size=9), legend.key.height=unit(0.15, "in"), legend.box.margin=unit(c(0.01,0.05,0.01,0.05), "in"), plot.margin=unit(c(0.05,0.1,0.1,0.1),"in"))


{cairo_pdf("output/figures/Fig03_spartials.pdf", width=6.5, height=9)

ggdraw() + draw_plot(goodobs, 0, 0.65, 0.5, 0.34)  + draw_plot(badobs, 0.5, 0.65, 0.5, 0.34) + draw_plot(spartials, 0.01, 0, 0.995, 0.65) + draw_plot_label(label=c("A", "B", "C", "D"), x=c(0, 0, 0.5, 0.48), y=c(0.995, 0.65, 0.995, 0.65))

}
dev.off()


#-------------------------------------------------------------------------
# Partials and spartials in each year with observations, for SI

sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range.shp")

states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_admin_1_states_provinces")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")



for(yr in 2008:2022){

# yr <- 2008

obs <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=states, fill="antiquewhite2", color="antiquewhite4") + 
	geom_sf(data=sdm.pres, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
	geom_point(data=filter(flow2, year==yr), aes(x=lon, y=lat, shape=flr, color=flr, size=flr), alpha=0.75) +
	scale_shape_manual(values=c(21,20), labels=c("No flowering", "Flowering"), name="Flowering") +
	scale_size_manual(values=c(1.1,1.3), labels=c("No flowering", "Flowering"), name="Flowering") +
	scale_color_manual(values=c("#ffffcc", "#253494"), labels=c("No flowering", "Flowering"), name="Flowering") + # nb these are false, true
	coord_sf(xlim = c(-119, -112.75), ylim = c(32.75, 38), expand = FALSE) +
	labs(title=paste(yr, "observations")) +
	theme_bw(base_size=11) + theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "slategray3"), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), legend.position=c(0.5, 0.1), legend.box.margin=unit(c(0.005, 0.005, 0.005, 0.005), "in"), legend.title=element_blank(), legend.text=element_text(size=9), plot.margin=unit(c(0.05,0.05,0.05,0.05), "in"), legend.key.size=unit(0.075, "in"), legend.key=element_rect(fill="#A8A378"), legend.spacing=unit(c(0.005, 0.005, 0.005, 0.005), "in"), legend.direction="horizontal", legend.background=element_rect(color="black", linewidth=0.1))

spart <- brick(paste("output/BART/Jotr_BART_spartials_", yr, ".grd", sep=""))
projection(spart)<-CRS("+init=epsg:4269")

spart.mask <- mask(spart, st_transform(sdm.pres[,2], crs=4269))

spart.df <- cbind(coordinates(spart.mask), as.data.frame(spart.mask)) %>% filter(!is.na(pptY1Y2)) |> rename(lon=x, lat=y) |> pivot_longer(all_of(jotr.preds), names_to="predictor", values_to="prFL") |> mutate(predictor=factor(predictor,jotr.preds))

levels(spart.df$predictor) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

spartplot <- ggplot(spart.df, aes(x=lon, y=lat, fill=prFL)) + geom_tile() + facet_wrap("predictor", nrow=3, labeller="label_parsed") + scale_fill_gradient(low="#ffffcc", high="#253494", name="Marginal Pr(Flowers)", breaks=c(0.5,0.55,0.6,0.65), limits=c(0.5,0.675), labels=c("", 0.55, 0.6, 0.65)) + labs(title=paste("Spatial partial effects for", yr)) + theme_minimal() + theme(axis.title=element_blank(), axis.text=element_blank(), legend.position="bottom", legend.text=element_text(size=9), panel.background=element_rect(fill="white"), panel.grid=element_blank())



{cairo_pdf(paste("output/figures/SIFigX_spartials_", yr, ".pdf", sep=""), width=4, height=9)

plot <- ggdraw() + draw_plot(obs, 0, 0.65, 1, 0.34) + draw_plot(spartplot, 0.05, 0, 0.85, 0.65) + draw_plot_label(label=c("A", "B"), x=0, y=c(0.995, 0.65))

print(plot)

}
dev.off()


} # END loop over years 




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

# LOOvalid <- read.csv("output/BART/year-year-LOO.csv")

mean(LOOvalid$AUC) # 0.60
sd(LOOvalid$AUC) # 0.18
sd(LOOvalid$AUC)/sqrt(nrow(LOOvalid)) # SE = 0.05

ggplot(LOOvalid, aes(x=N_flr, y=AUC)) + geom_point()
ggplot(LOOvalid, aes(x=N_noflr, y=AUC)) + geom_point()
ggplot(LOOvalid, aes(x=N_flr+N_noflr, y=AUC)) + geom_point()
ggplot(LOOvalid, aes(x=N_flr/N_noflr, y=AUC)) + geom_point()




#-------------------------------------------------------------------------
# "recent" model based only on years with >100 observations

table(flow2$year) # at least 100/yr from 2016 onward

flow4 <- flow2 |> filter(year >= 2016)

jotr.recent.mod <- bart(y.train=as.numeric(flow4[,"flr"]), x.train=flow4[,jotr.preds], keeptrees=TRUE)

summary(jotr.recent.mod)

invisible(jotr.recent.mod$fit$state)
write_rds(jotr.recent.mod, file="output/BART/bart.recent.model.Jotr.rds")
# jotr.mod <- read_rds("output/BART/bart.recent.model.Jotr.rds")

# cross-validate with earlier data ---------------
earlyValid <- data.frame(matrix(0,0,5))
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

auc # 0.65, okay


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

# LOOvalidRI <- read.csv( "output/BART/RI_year-year-LOO.csv")

mean(LOOvalidRI$AUC) # 0.59
sd(LOOvalidRI$AUC) # 0.15
sd(LOOvalidRI$AUC)/sqrt(nrow(LOOvalidRI)) # SE = 0.04

#-------------------------------------------------------------------------
# models fit to subspecies data

# eastern first
yuja.mod <- bart(y.train=as.numeric(yuja[,"flr"]), x.train=yuja[,jotr.preds], keeptrees=TRUE)

summary(yuja.mod)

invisible(yuja.mod$fit$state)
write_rds(yuja.mod, file="output/BART/bart.model.yuja.rds")
# jotr.mod <- read_rds("output/BART/bart.model.yuja.rds")

# then western
yubr.mod <- bart(y.train=as.numeric(yubr[,"flr"]), x.train=yubr[,jotr.preds], keeptrees=TRUE)

summary(yubr.mod)

invisible(yubr.mod$fit$state)
write_rds(yubr.mod, file="output/BART/bart.model.yubr.rds")
# jotr.mod <- read_rds("output/BART/bart.model.yubr.rds")


