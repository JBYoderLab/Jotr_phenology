# Comparing predicted historical flowering to independent validation records
# Assumes local environment
# jby 2024.02.14

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("raster")
library("sp")
library("sf")
library("hexbin")
library("embarcadero")

source("../shared/Rscripts/base.R") # my special mix of personal functions
source("../shared/Rscripts/base_graphics.R") # my special mix of personal functions

library("ggdark")

#-------------------------------------------------------------------------
# load up and organize historical data

# Jotr predicted flowering

jotr.histStack <- raster::stack("output/BART/jotr_BART_predicted_flowering_1900-2023_nomask.grd")
jotr.ri.histStack <- raster::stack("output/BART/jotr_BART_RI_predicted_flowering_1900-2023_nomask.grd")

yubr.histStack <- raster::stack("output/BART/YUBR_BART_predicted_flowering_1900-2023_nomask.grd")
yuja.histStack <- raster::stack("output/BART/YUJA_BART_predicted_flowering_1900-2023_nomask.grd")


# expert validation observations
vobs <- read.csv("data/Validation_obs_by_spp.csv")

glimpse(vobs)

table(vobs$obs_by)


#-------------------------------------------------------------------------
# Incorporating NPN records

# US NPN data
npn <- read.csv("data/NPN_datasheet_1704758997752/status_intensity_observation_data.csv")

glimpse(npn)

table(npn$Partner_Group) # JTNP is the single external source
length(unique(npn$Site_ID)) # 16 unique sites
length(unique(paste(npn$Longitude, npn$Latitude))) # confirm 16 locations at resolution of this data

table(year(ymd(npn$Observation_Date)), npn$Site_ID) # records for 13 years, 2011-2023
sum(table(year(ymd(npn$Observation_Date)), npn$Site_ID)>0) # works out to 108 location-year combos

table(npn$Phenophase_Category) # ONLY flowers or fruits recorded

table(npn$Intensity_Category_ID)
# 23 = fruits present
# 24 = ripe fruits present
# 25 = recent fruit drop
# 31 = open flowers (peak)
# 35 = flowers or flower heads present
# 48 = flowers and flower buds 
# 50 = open flowers percentage (individual)
# 56 = fruits present
# 58 = ripe fruit percentage
# 59 = recent fruit or seed drop

# How are these organized?
filter(npn, year(ymd(Observation_Date))==2011, Site_ID==5389) %>% dplyr::select(Observation_Date, Phenophase_Category)

# where are these located?
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range.shp")
ggplot() + geom_sf(data=sdm.pres) + geom_point(data=npn, aes(x=Longitude, y=Latitude))

filter(npn, Longitude < -122) # HUH
filter(npn, Longitude > -113) %>% group_by(year(ymd(Observation_Date))) %>% summarize(n=length(Observation_Date))
filter(npn, Longitude > -113, year(ymd(Observation_Date))==2018) %>% group_by(Phenophase_Category) %>% summarize(n=length(Observation_Date))

# okay, let's try this ... natural range records only
npn.vobs <- npn %>% filter(Longitude > -120, Longitude < -114) %>% 
				   mutate(year = year(ymd(Observation_Date))) %>% 
				   filter(!duplicated(paste(Site_ID, year))) %>%
				   dplyr::select(Longitude, Latitude, Site_Name, year) %>%
				   rename(lon=Longitude, lat=Latitude, location=Site_Name) %>%
				   mutate(type="YUBR", obs_by="NPN", obs_flowers=TRUE, obs_fruit=TRUE, obs_moths=NA, obs_no_flowers=NA)

glimpse(npn.vobs)


# for comparison to our data ...
flow2 <- read.csv("output/flowering_obs_climate_subsp.csv") %>% filter(!(year==2018.5 & flr==TRUE), year>=2008) %>% mutate(year=floor(year)) # drop the late-flowering anomaly

glimpse(flow2) # 2,632 in our final working set
length(unique(paste(flow2$lat, flow2$lon))) # 1140

#-------------------------------------------------------------------------
# Incorporating annotations from MDEP survey plot photos

mdep <- read.csv("output/MDEP_survey_records_cleaned.csv") %>% mutate(obs_flowers=flr, obs_fruit=flr)

glimpse(mdep)

#-------------------------------------------------------------------------
# Incorporating pre-2008 iNat records

inatEarly <- read.csv("output/flowering_obs_climate_subsp.csv") %>% filter(!(year==2018.5 & flr==TRUE), year<2008) %>% mutate(year=floor(year)) %>% dplyr::select(lat, lon, year, flr, type) %>% mutate(obs_by="iNat", location=NA, obs_flowers=NA, obs_fruit=NA, obs_moths=NA, obs_no_flowers=NA)


#-------------------------------------------------------------------------
# Incorporating CCH records

herb <- read.csv("data/SymbOutput_2024-02-11_233150_DwC-A/occurrences_cleaned.csv")

herb$type <- NA
herb$type[grepl(".*brevifolia.*", herb$scientificName)] <- "YUBR"
herb$type[grepl(".*jaegeriana.*", herb$scientificName)] <- "YUJA"

glimpse(herb)
table(herb$year) # LMAO back to 1876??
table(herb$type)
hist(herb$coordinateUncertaintyInMeters) # hoo boy

herb.vobs <- herb %>% filter(year>=1900, coordinateUncertaintyInMeters<=4000) %>% 
	mutate(obs_flowers=NA, obs_fruit=NA, obs_moths=NA, obs_no_flowers=NA) %>% 
	dplyr::select(lat, lon, location, type, year, obs_by, obs_flowers, obs_fruit, obs_moths, obs_no_flowers, flr)

glimpse(herb.vobs)

#-------------------------------------------------------------------------
# Merge validation sources

vobs.all <- rbind(vobs, npn.vobs) %>% mutate(flr=(obs_flowers | obs_fruit)) %>% rbind(herb.vobs) #%>% rbind(inatEarly) # %>% rbind(inat)

glimpse(vobs.all)
table(vobs.all$obs_by, vobs.all$year)
table(vobs.all$flr, vobs.all$year)
table(vobs.all$flr, vobs.all$obs_by)

#-------------------------------------------------------------------------
# Link validation records to predictions

# best-power prediction thresholds from the original models
# jotr: 0.26
# jotr RI: 0.25

# YUBR: 0.28
# YUJA: 0.16

glimpse(vobs.all) 

valid <- NULL # initialize, in a lazy way

for(yr in sort(unique(vobs.all$year))){

# yr <- 2023

vsub <- filter(vobs.all, year==yr)

vsub$jotr_prFlr <- raster::extract(jotr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$jotr_Flr <- vsub$jotr_prFlr >= 0.26

vsub$jotr_RI.prFlr <- raster::extract(jotr.ri.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$jotr_RI.Flr <- vsub$jotr_RI.prFlr >= 0.25


vsub$YUBR_prFlr <- raster::extract(yubr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$YUBR_Flr <- vsub$YUBR_prFlr >= 0.28

vsub$YUJA_prFlr <- raster::extract(yuja.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$YUJA_Flr <- vsub$YUJA_prFlr >= 0.16

valid <- rbind(vsub, valid)

}

valid <- filter(valid, !is.na(jotr_prFlr), year < 2023)

glimpse(valid) 
table(valid$obs_by, useNA="ifany")
head(valid) # cool
tail(valid)

write.table(valid, "output/expert_obs_validations.csv", sep=",", col.names=TRUE, row.names=FALSE)

# valid <- read.csv("output/expert_obs_validations.csv")

# RANGE-WIDE -------------

# without RI
performance(prediction(valid$jotr_prFlr, valid$flr), "auc")@y.values[[1]] # 0.69 overall

t.test(jotr_prFlr~flr, data=valid) # p = 2.885e-05 okay!

performance(prediction(filter(valid, obs_by=="CCH")$jotr_prFlr, filter(valid, obs_by=="CCH")$flr), "auc")@y.values[[1]]
performance(prediction(filter(valid, obs_by!="CCH")$jotr_prFlr, filter(valid, obs_by!="CCH")$flr), "auc")@y.values[[1]]


# jotr, with RI
performance(prediction(valid$jotr_RI.prFlr, valid$flr), "auc")@y.values[[1]] # 0.66 overall

t.test(jotr_RI.prFlr~flr, data=valid) # p = 7.567e-07

# SINGLE SPECIES -------------
vYUBR <- filter(valid, type=="YUBR")
vYUJA <- filter(valid, type=="YUJA")

# YUBR
performance(prediction(vYUBR$YUBR_prFlr, vYUBR$flr), "auc")@y.values[[1]] # 0.65 overall, whee

t.test(YUBR_prFlr~flr, data=vYUBR) # p = 0.0003059, okay

# YUJA
performance(prediction(vYUJA$YUJA_prFlr, vYUJA$flr), "auc")@y.values[[1]] # 0.56 overall, whee

t.test(YUJA_prFlr~flr, data=vYUJA) # p = 0.06305, fine


#-------------------------------------------------------------------------
# newspaper records

news <- read.csv("data/news_accounts/news_reports.csv")
glimpse(news)
table(news$year)

# For a first attempt, let's just look at range-wide average ...

jotr.flr <- read.csv("output/historic_flowering_reconst_jotr.csv")
glimpse(jotr.flr)

# What we really need, though, is to use what geographic information we have
jotr.histStack <- raster::stack("output/BART/jotr_BART_predicted_flowering_1900-2023.grd")

# relevant polygons
JTNP <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_parks_and_protected_lands_scale_rank") %>% filter(unit_code=="JOTR") # the national park, when that's specified
Lanc <- read_sf(dsn = "../data/spatial/CA_Counties/", lay= "CA_Counties_TIGER2016") %>% filter(NAME=="Los Angeles") # LA County as a proxy for "Lancaster" or "Cajon" or similar
Kern <- read_sf(dsn = "../data/spatial/CA_Counties/", lay= "CA_Counties_TIGER2016") %>% filter(NAME=="Kern") # and Kern

# then mask
JTNP.histFlr <- mask(jotr.histStack, st_transform(JTNP, crs=4269), touches=TRUE) # go generous
Lanc.histFlr <- mask(jotr.histStack, st_transform(Lanc, crs=4269), touches=TRUE) # go generous
Kern.histFlr <- mask(jotr.histStack, st_transform(Kern, crs=4269), touches=TRUE) # go generous

# then this gets complex
news.prFLs <- news %>% mutate(prFL.loc=NA, prFL.all=NA)

for(i in 1:nrow(news.prFLs)){

if(news.prFLs$location[i]=="JTNP") loc.yr <- JTNP.histFlr[[paste("prFL", news.prFLs$year[i], sep=".")]]
if(news.prFLs$location[i]=="Lancaster") loc.yr <- Lanc.histFlr[[paste("prFL", news.prFLs$year[i], sep=".")]]
if(news.prFLs$location[i]=="Kern County") loc.yr <- Kern.histFlr[[paste("prFL", news.prFLs$year[i], sep=".")]]

news.prFLs$prFL.all[i] <- cellStats(jotr.histStack[[paste("prFL", news.prFLs$year[i], sep=".")]], "mean")
news.prFLs$prFL.loc[i] <- cellStats(loc.yr, "mean")

}

glimpse(news.prFLs)

write.table(news.prFLs, "output/news_reports_prFL.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

t.test(prFL.loc~flr, data=news.prFLs) # p = 0.0004, YES
t.test(prFL.all~flr, data=news.prFLs) # p = 0.003


news.prFLs %>% group_by(flr) %>% summarize(mnPrLoc = mean(prFL.loc), sePrLoc = sd(prFL.loc)/sqrt(length(prFL.loc)), mnPrAll = mean(prFL.all), sePrAll = sd(prFL.all)/sqrt(length(prFL.all)))



ggplot() + geom_histogram(data=prFL.sum, aes(x=mn.prFL), fill="gray60") + geom_vline(data=news.prFLs, aes(xintercept=prFL.loc, color=flr)) + geom_vline(xintercept=0.26, linetype=2)


# visualize news reports and other validation sources ................
valid_srcs <- valid |> dplyr::select(obs_by, jotr_prFlr, flr) |> rename(fitted=jotr_prFlr, observed=flr) |> mutate(observed=as.numeric(observed), type="Validation records", classified=fitted>0.2635843)
glimpse(valid_srcs)

accu <- rbind(obs_valid, exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.2635843)

ggplot(valid_srcs, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.26, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + 

facet_wrap("obs_by", ncol=1) + 

annotate(geom="text", x=0.27, y=1.45, label="Classfication cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))


{cairo_pdf("output/figures/FigXX_validation_sources.pdf", width=8, height=6)

ggdraw() + draw_plot(mod_varimp, 0, 0.53, 0.45, 0.44) + draw_plot(mod_accuracy, 0, 0, 0.45, 0.55) + draw_plot(partplot, 0.45, 0, 0.55, 1) + draw_plot_label(label=c("A", "B", "C"), x=c(0, 0, 0.45), y=c(1,0.55,1))

}
dev.off()



#-------------------------------------------------------------------------
# visualize accuracy with original data and validations 

jotr.mod <- read_rds("output/BART/bart.model.Jotr.rds")
names(jotr.mod)

obs_valid <- summary(jotr.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data")
glimpse(obs_valid) 

exp_valid <- valid |> dplyr::select(jotr_prFlr, flr) |> rename(fitted=jotr_prFlr, observed=flr) |> mutate(observed=as.numeric(observed), type="Validation records")
glimpse(exp_valid)

accu <- rbind(obs_valid, exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.2635843)

# model accuracy figure .........................
mod_accuracy <- ggplot(accu, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.26, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + facet_wrap("type", nrow=2) + 

annotate(geom="text", x=0.27, y=1.45, label="Classfication cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

# and then also the predictor selection figure ............
jotr.varimp <- read_rds("output/BART/bart.varimp.Jotr.rds")

jotr.varimp$data <- jotr.varimp$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

levels(jotr.varimp$data$names) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

jotr.varimp$labels$group <- "Trees"
jotr.varimp$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

mod_varimp <- ggplot(jotr.varimp$data, aes(x=names, y=imp, color=trees, group=trees)) + geom_line(linewidth=0.75) + geom_point(size=1.5, color="gray30") + scale_color_manual(values=c('#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84'), name="Trees") + labs(y="Relative contribution") + scale_x_discrete(label=label_parse) + theme_minimal() + theme(legend.position=c(0.8, 0.75), axis.text=element_text(size=9), axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank(), legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size=unit(0.15, "in"))  # okay nice

# and finally the partial effects ...............
p <- read_rds("output/BART/bart.model.Jotr.partials.rds")

partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

partplot <- ggplot(partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

partplot

# put it all togeether
{cairo_pdf("output/figures/Fig02_predictors_accuracy_partials.pdf", width=8, height=6)

ggdraw() + draw_plot(mod_varimp, 0, 0.53, 0.45, 0.44) + draw_plot(mod_accuracy, 0, 0, 0.45, 0.55) + draw_plot(partplot, 0.45, 0, 0.55, 1) + draw_plot_label(label=c("A", "B", "C"), x=c(0, 0, 0.45), y=c(1,0.55,1))

}
dev.off()


# figure to do the above but with the two species separately -------------

jotr.preds <- jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

# first YUJA 
yuja.mod <- read_rds("output/BART/bart.model.yuja.rds")
names(yuja.mod)

yuja_obs_valid <- summary(yuja.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data")
glimpse(yuja_obs_valid) 

yuja_exp_valid <- valid |> filter(type=="YUJA", !is.na(obs_flowers)) |> dplyr::select(YUJA_prFlr, obs_flowers) |> rename(fitted=YUJA_prFlr, observed=obs_flowers) |> mutate(observed=as.numeric(observed), type="Validation records")
glimpse(yuja_exp_valid)

yuja_accu <- rbind(yuja_obs_valid, yuja_exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.1589413)

# model accuracy figure
yuja_accuracy <- ggplot(yuja_accu, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.16, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + facet_wrap("type", nrow=2) + 

annotate(geom="text", x=0.17, y=1.45, label="Classfication cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

labs(title="Model accuracy, YUJA") + theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

yuja_accuracy

yuja_p <- partial(yuja.mod, jotr.preds, trace=FALSE, smooth=5)

yuja_partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", yuja_p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", yuja_p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", yuja_p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", yuja_p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", yuja_p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", yuja_p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

yuja_partplot <- ggplot(yuja_partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)", title="Partial effects, YUJA") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

yuja_partplot


# now YUBR
yubr.mod <- read_rds("output/BART/bart.model.yubr.rds")
names(yubr.mod)

yubr_obs_valid <- summary(yubr.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data")
glimpse(yubr_obs_valid) 

yubr_exp_valid <- valid |> filter(type=="YUBR", !is.na(obs_flowers)) |> dplyr::select(YUBR_prFlr, obs_flowers) |> rename(fitted=YUBR_prFlr, observed=obs_flowers) |> mutate(observed=as.numeric(observed), type="Validation records")
glimpse(yubr_exp_valid)

yubr_accu <- rbind(yubr_obs_valid, yubr_exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.2776204)

# model accuracy figure
yubr_accuracy <- ggplot(yubr_accu, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.28, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + facet_wrap("type", nrow=2) + 

annotate(geom="text", x=0.29, y=1.45, label="Classfication cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

labs(title="Model accuracy, YUBR") + theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

yubr_accuracy

yubr_p <- partial(yubr.mod, jotr.preds, trace=FALSE, smooth=5)

yubr_partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", yubr_p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", yubr_p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", yubr_p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", yubr_p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", yubr_p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", yubr_p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

yubr_partplot <- ggplot(yubr_partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)", title="Partial effects, YUBR") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

yubr_partplot


{cairo_pdf("output/figures/SIFig_subspecies_accuracy_partials.pdf", width=6.5, height=8)

ggdraw() + draw_plot(yubr_accuracy, 0, 0.6, 0.5, 0.4) + draw_plot(yuja_accuracy, 0.5, 0.6, 0.5, 0.4) + draw_plot(yubr_partplot, 0, 0, 0.5, 0.6) + draw_plot(yuja_partplot, 0.5, 0, 0.5, 0.6) + draw_plot_label(label=c("A", "B", "C", "D"), x=c(0, 0, 0.5, 0.5), y=c(1, 0.6, 1, 0.6))

}
dev.off()

