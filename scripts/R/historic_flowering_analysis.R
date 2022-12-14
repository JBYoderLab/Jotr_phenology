# Analyzing predicted historical flowering in Joshua tree
# Assumes local environment
# jby 2022.11.22

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("raster")
library("sp")
library("sf")
library("hexbin")


source("../shared/Rscripts/base.R") # my special mix of personal functions
source("../shared/Rscripts/base_graphics.R") # my special mix of personal functions

library("ggdark")

#-------------------------------------------------------------------------
# load up and organize historical data

# Jotr SDM and species boundaries
sdm.pres <- read_sf("../data/Yucca/jotr_BART_sdm_pres", "jotr_BART_sdm_pres")
spp.ranges <- read_sf("data/Jotr_ssp_range.kml")

yubr.pres <- st_zm(st_intersection(sdm.pres, st_transform(filter(spp.ranges, Name=="YUBR"), crs=crs(sdm.pres))))
yuja.pres <- st_zm(st_intersection(sdm.pres, st_transform(filter(spp.ranges, Name=="YUJA"), crs=crs(sdm.pres))))

# flowering observation data
obs <- read.csv("output/flowering_obs_climate_v2_subsp.csv") %>% mutate(y2 = year) %>% mutate(year=floor(year), type = factor(type, c("YUBR", "YUJA")))

# raster files of predicted prFL
jotr.files <- list.files("output/BART/predictions", pattern=".bil", full=TRUE)
yubr.files <- list.files("output/BART/predictions.YUBR", pattern=".bil", full=TRUE)
yuja.files <- list.files("output/BART/predictions.YUJA", pattern=".bil", full=TRUE)

MojExt <- extent(-119, -112, 33, 38) # Mojave extent, maybe useful

# expert validation observations
vobs <- read.csv("data/Validation_obs_by_spp.csv")


# Jotr --------------------------------
jotr.histStack <- raster::stack(sapply(jotr.files, function(x) crop(raster::raster(x), MojExt)))
names(jotr.histStack) <- paste("prFL",1900:2022,sep=".")
projection(jotr.histStack)<-CRS("+init=epsg:4269")

jotr.histStack

jotr.maskHist <- mask(jotr.histStack, st_transform(sdm.pres[,2], crs=4269))

jotr.hist.flowering <- cbind(coordinates(jotr.maskHist), as.data.frame(jotr.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(jotr.hist.flowering)

write.table(jotr.hist.flowering, "output/historic_flowering_reconst_jotr.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)



# YUBR --------------------------------
yubr.histStack <- raster::stack(sapply(yubr.files, function(x) crop(raster::raster(x), MojExt)))
names(yubr.histStack) <- paste("prFL",1900:2022,sep=".")
projection(yubr.histStack)<-CRS("+init=epsg:4269")

yubr.histStack

yubr.maskHist <- mask(yubr.histStack, st_transform(yubr.pres[,2], crs=4269))

yubr.hist.flowering <- cbind(coordinates(yubr.maskHist), as.data.frame(yubr.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(yubr.hist.flowering)

write.table(yubr.hist.flowering, "output/historic_flowering_reconst_YUBR.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


# YUJA --------------------------------
yuja.histStack <- raster::stack(sapply(yuja.files, function(x) crop(raster::raster(x), MojExt)))
names(yuja.histStack) <- paste("prFL",1900:2022,sep=".")
projection(yuja.histStack)<-CRS("+init=epsg:4269")

yuja.histStack

yuja.maskHist <- mask(yuja.histStack, st_transform(yuja.pres[,2], crs=4269))

yuja.hist.flowering <- cbind(coordinates(yuja.maskHist), as.data.frame(yuja.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(yuja.hist.flowering)

write.table(yuja.hist.flowering, "output/historic_flowering_reconst_YUJA.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


# YUJA+YUBR modeled separately --------------------------------
hist.flowering <- rbind(data.frame(type="YUBR", yubr.hist.flowering), data.frame(type="YUJA", yuja.hist.flowering))

write.table(hist.flowering, "output/historic_flowering_reconst_YUBR-YUJA.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

#-------------------------------------------------------------------------
# Validation observations

# best-power prediction thresholds from the original models
# jotr: 0.26
# YUBR: 0.30
# YUJA: 0.23

glimpse(vobs)

valid <- NULL # initialize, in a lazy way

for(yr in unique(vobs$year)){

# yr <- 2005

vsub <- filter(vobs, year==yr)

vsub$jotr_prFlr <- raster::extract(jotr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$jotr_Flr <- raster::extract(jotr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.26

vsub$YUBR_prFlr <- raster::extract(yubr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$YUBR_Flr <- raster::extract(yubr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.30

vsub$YUJA_prFlr <- raster::extract(yuja.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$YUJA_Flr <- raster::extract(yuja.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.23

valid <- rbind(vsub, valid)

}

head(valid) # cool

write.table(valid, "output/expert_obs_validations.csv", sep=",", col.names=TRUE, row.names=FALSE)


# species-specific models/validation

table(valid$obs_flowers, valid$jotr_Flr) # eesh
chisq.test(table(valid$obs_flowers, valid$jotr_Flr)) # LOLsob

t.test(jotr_prFlr~obs_flowers, data=valid) # n.s., dammit

vYUBR <- filter(valid, type=="YUBR")
vYUJA <- filter(valid, type=="YUJA")

t.test(YUBR_prFlr~obs_flowers, data=vYUBR) # p = 0.051 LOL
t.test(YUJA_prFlr~obs_flowers, data=vYUJA) # n.s.
chisq.test(table(vYUBR$obs_flowers, vYUBR$YUBR_Flr)) # LOLsob
chisq.test(table(vYUJA$obs_flowers, vYUJA$YUJA_Flr)) # LOLsob

table(c(vYUBR$obs_flowers,vYUJA$obs_flowers), c(vYUBR$YUBR_Flr, vYUJA$YUJA_Flr)) # eesh
chisq.test(table(c(vYUBR$obs_flowers,vYUJA$obs_flowers), c(vYUBR$YUBR_Flr, vYUJA$YUJA_Flr))) # LOLsob

# just RY observations of YUBR

t.test(YUBR_prFlr~obs_flowers, data=filter(vYUBR, obs_by=="Ray Yeager"))
t.test(YUBR_prFlr~obs_flowers, data=filter(vYUBR, obs_by=="CIS")) # p = 0.007053

chisq.test(table(filter(vYUBR, obs_by=="CIS")$obs_flowers, filter(vYUBR, obs_by=="CIS")$YUBR_Flr)) # HMM


#-------------------------------------------------------------------------
# observations versus same-year predictions

hist.flowering <- read.csv("output/historic_flowering_reconst_YUBR-YUJA.csv", h=TRUE) %>% mutate(type = factor(type, c("YUBR", "YUJA")))

glimpse(hist.flowering)

{cairo_pdf("output/figures/obs-vs-prediction_2021.pdf", width=9, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==2021), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.5)) + 

geom_point(data=filter(obs, year==2021), aes(x=lon, y=lat, shape=flr, color=flr), size=1) + 
scale_shape_manual(values=c(1,20), name="Flowers seen?") + 
scale_color_manual(values=c("gray60", park_palette("JoshuaTree")[2]), name="Flowers seen?") +

labs(x="Latitude", y="Longitude") +

facet_grid(.~type, scale="free_x") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.9,0.1,0.9), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=8), plot.background=element_rect(color="black", fill="black"))

}
dev.off()

{cairo_pdf("output/figures/obs-vs-prediction_2022.pdf", width=9, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==2022), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.5)) + 

geom_point(data=filter(obs, year==2022), aes(x=lon, y=lat, shape=flr, color=flr), size=1) + 
scale_shape_manual(values=c(1,20), name="Flowers seen?") + 
scale_color_manual(values=c("gray60", park_palette("JoshuaTree")[2]), name="Flowers seen?") +

labs(x="Latitude", y="Longitude") +

facet_grid(.~type, scale="free_x") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.9,0.1,0.9), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=8), plot.background=element_rect(color="black", fill="black"))

}
dev.off()

# all years with observations, YUBR
{cairo_pdf("output/figures/obs-vs-prediction_YUBR.pdf", width=9.5, height=7.5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year%in%2010:2022, type=="Western"), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.5)) + 

geom_point(data=filter(obs, type=="Western"), aes(x=lon, y=lat, shape=flr, color=flr), size=0.75) + 
scale_shape_manual(values=c(1,20), name="Flowers seen?") + 
scale_color_manual(values=c("gray60", park_palette("JoshuaTree")[2]), name="Flowers seen?") +

facet_wrap("year", nrow=3) + 

dark_mode(theme_minimal()) + theme(legend.position=c(0.85,0.2), legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.title=element_blank(), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()

# all years with observations, YUJA
{cairo_pdf("output/figures/obs-vs-prediction_YUJA.pdf", width=9.5, height=7.5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year%in%2010:2022, type=="Eastern"), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.5)) + 

geom_point(data=filter(obs, type=="Eastern"), aes(x=lon, y=lat, shape=flr, color=flr), size=0.75) + 
scale_shape_manual(values=c(1,20), name="Flowers seen?") + 
scale_color_manual(values=c("gray60", park_palette("JoshuaTree")[2]), name="Flowers seen?") +

facet_wrap("year", nrow=3) + 

dark_mode(theme_minimal()) + theme(legend.position=c(0.85,0.2), legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.title=element_blank(), axis.text=element_blank(), plot.margin=unit(c(0.5,0.1,0.75,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()



#-------------------------------------------------------------------------
# Individual years

# my birth year, LOL
{cairo_pdf("output/figures/prediction_1982.pdf", width=5.5, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==1982), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()

# the year of McKelvey's description of YUJA
{cairo_pdf("output/figures/prediction_1935.pdf", width=5.5, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==1935), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


# year Joshua tree National Monument established
{cairo_pdf("output/figures/prediction_1936.pdf", width=5.5, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==1936), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


#-------------------------------------------------------------------------
# TRENDS within cells


# cor(flr, yr) --- crude but it's a start -----------------

cellcors <- jotr.hist.flowering %>% group_by(lon,lat) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman"))) # rangewide model

cellcors.sep <- hist.flowering %>% group_by(lon,lat,type) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman"))) # two types together

glimpse(cellcors)
glimpse(cellcors.sep)


{cairo_pdf("output/figures/prFL-vs-time_correlations_YUBR-YUJA.pdf", width=3, height=4)

ggplot(cellcors.sep, aes(x=estimate, fill=p.value<=0.01, color=p.value<=0.01)) + geom_histogram(color=NA) + 

scale_fill_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +
scale_color_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +

#annotate("text", x=0.1, y=150, label="All cells", color=park_palette("JoshuaTree")[7], fontface="bold", family="Arial Narrow") + 

#annotate("text", x=-0.15, y=50, label="Cells with\ncorrelation\np < 0.01", color=park_palette("JoshuaTree")[2], fontface="bold", family="Arial Narrow", lineheight=0.8) + 

geom_vline(xintercept=0) +

facet_grid(type~., scale="free_x") +

labs(x=expression("Spearman's"~rho), y="Grid cells") + dark_mode(theme_minimal()) + theme(legend.position="none", plot.background=element_rect(color="black", fill="black"), plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()

{cairo_pdf("output/figures/prFL-vs-time_correlations_jotr.pdf", width=3, height=2.5)

ggplot(cellcors, aes(x=estimate, fill=p.value<=0.01, color=p.value<=0.01)) + geom_histogram(color=NA) + 

scale_fill_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +
scale_color_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +

#annotate("text", x=0.1, y=150, label="All cells", color=park_palette("JoshuaTree")[7], fontface="bold", family="Arial Narrow") + 

#annotate("text", x=-0.15, y=50, label="Cells with\ncorrelation\np < 0.01", color=park_palette("JoshuaTree")[2], fontface="bold", family="Arial Narrow", lineheight=0.8) + 

geom_vline(xintercept=0) +

labs(x=expression("Spearman's"~rho), y="Grid cells") + dark_mode(theme_minimal()) + theme(legend.position="none", plot.background=element_rect(color="black", fill="black"), plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()

# Plot it on a map
{cairo_pdf("output/figures/prFL-vs-time_map_jotr.pdf", width=5.5, height=6)

ggplot(cellcors, aes(x=lon, y=lat, fill=estimate)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="div", palette=1, direction=1, name=expression("Spearman's"~rho)) + labs(x="Longitude", y="Latitude", title="Correlation between probability of flowering and time") + 

dark_mode(theme_minimal(base_size=12)) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black"))

}
dev.off()

{cairo_pdf("output/figures/prFL-vs-time_map_YUBR-YUJA.pdf", width=5.5, height=6)

ggplot(cellcors.sep, aes(x=lon, y=lat, fill=estimate)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="div", palette=1, direction=1, name=expression("Spearman's"~rho)) + labs(x="Longitude", y="Latitude", title="Correlation between probability of flowering and time") + 

dark_mode(theme_minimal(base_size=12)) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black"))

}
dev.off()


# number of flowering years/decade? -----------------------

# best-power prediction thresholds from the original models
# jotr: 0.26
# YUBR: 0.30
# YUJA: 0.23

flyrs.jotr <- jotr.hist.flowering %>% group_by(lon,lat) %>% summarize(flyrs=length(which(prFL>=0.26))) # rangewide model
range(jotr.hist.flowering$year) # remember: 1900-2022
glimpse(flyrs.jotr)

{cairo_pdf("output/figures/flowering-years_jotr.pdf", width=3.5, height=3.5)
ggplot(flyrs.jotr, aes(x=flyrs)) + geom_histogram() + labs(x = "Years with predicted flowering, 1900-2022", y="Cells")
}
dev.off()



ggplot(jotr.hist.flowering, aes(x=prFL)) + geom_histogram() + labs(x = "Modeled pr(FL)", y="Cells")


flyrs.yubr <- hist.flowering |> filter(type=="YUBR") |> group_by(lon,lat) |> summarize(flyrs=length(which(prFL>=0.30))) # YUBR
flyrs.yuja <- hist.flowering |> filter(type=="YUJA") |> group_by(lon,lat) |> summarize(flyrs=length(which(prFL>=0.23))) # YUJA

flyrs.sep <- rbind(data.frame(type="YUBR", flyrs.yubr), data.frame(type="YUJA", flyrs.yuja))

{cairo_pdf("output/figures/flowering-years_YUBR-YUJA.pdf", width=6, height=3.5)
ggplot(flyrs.sep, aes(x=flyrs)) + geom_histogram() + facet_wrap("type") + labs(x = "Years with predicted flowering, 1900-2022", y="Cells")
}
dev.off()


ggplot(hist.flowering, aes(x=prFL)) + geom_histogram() + facet_wrap("type") + labs(x = "Modeled pr(flowering)", y="Cells")

# pre decade ----------------------------------------------

flyrs.dec.jotr <- jotr.hist.flowering %>% mutate(decade=paste(floor(year/10)*10,"-9", sep="")) %>% group_by(lon,lat,decade) %>% summarize(flyrs=length(which(prFL>=0.26)))

glimpse(flyrs.dec.jotr)

ggplot(flyrs.dec.jotr, aes(x=decade, y=flyrs)) + geom_boxplot() + labs(x="Decade", y="Predicted flowering years") # hrmmmm


flyrs.dec.yubr <- hist.flowering %>% filter(type=="YUBR") %>% mutate(decade=paste(floor(year/10)*10,"-9", sep="")) %>% group_by(lon,lat,decade) %>% summarize(flyrs=length(which(prFL>=0.3))) %>% mutate(type="YUBR")

flyrs.dec.yuja <- hist.flowering %>% filter(type=="YUJA") %>% mutate(decade=paste(floor(year/10)*10,"-9", sep="")) %>% group_by(lon,lat,decade) %>% summarize(flyrs=length(which(prFL>=0.23))) %>% mutate(type="YUJA")

flyrs.dec.sep <- rbind(flyrs.dec.yubr, flyrs.dec.yuja)

glimpse(flyrs.dec.sep)

ggplot(flyrs.dec.sep, aes(x=decade, y=flyrs, fill=type)) + geom_boxplot() + labs(x="Decade", y="Predicted flowering years") # HRMMMMM



# Plot it on a map
{cairo_pdf("output/figures/flyrs_map_jotr.pdf", width=5.5, height=6)

ggplot(flyrs.jotr, aes(x=lon, y=lat, fill=flyrs)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="div", palette=1, direction=1, name="Years flowering predicted") + labs(x="Longitude", y="Latitude", title="Predicted flowering years, 1900-2022") + 

dark_mode(theme_minimal(base_size=12)) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black"))

}
dev.off()

{cairo_pdf("output/figures/flyrs_map_YUBR-YUJA.pdf", width=5.5, height=6)

ggplot(flyrs.sep, aes(x=lon, y=lat, fill=flyrs)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="div", palette=1, direction=1, name="Years flowering predicted") + labs(x="Longitude", y="Latitude", title="Predicted flowering years, 1900-2022") + 

dark_mode(theme_minimal(base_size=12)) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black"))

}
dev.off()


#-------------------------------------------------------------------------
# rangewide trends

# make this usable ...
sub.hist.flowering <- hist.flowering %>% group_by(year) %>% slice_sample(prop=0.05) 

# with points
{png("output/figures/prFL-vs-time.png", width=1400, height=800)

ggplot(sub.hist.flowering, aes(x=year, y=prFL)) + geom_jitter(alpha=0.2, color=park_palette("JoshuaTree")[7]) + geom_smooth(method="loess", color=park_palette("JoshuaTree")[7]) + labs(x="Year", y="Flowering probability") + dark_mode(theme_ipsum(base_size=30, axis_title_size=36)) + theme(plot.background=element_rect(color="black")) # weeee

}
dev.off()

# try this ...

{cairo_pdf("output/figures/prFL-vs-time_hex.pdf", width=5, height=4)

ggplot(sub.hist.flowering, aes(x=year, y=prFL)) + geom_hex() + geom_hline(yintercept=median(filter(sub.hist.flowering, year<1930)$prFL), linetype=2) + geom_smooth(method="loess", span=0.75, color=park_palette("JoshuaTree")[5], se=FALSE) + annotate("text", x=2021, y=1.05*median(filter(sub.hist.flowering, year<1930)$prFL), label="Median, 1900-1929", hjust=1, vjust=0, size=3) + scale_fill_gradient(low="gray95", high=park_palette("JoshuaTree")[5]) + labs(x="Year", y="Probability of flowering") + theme_ipsum() + theme(legend.position="none")

}
dev.off()


# sumarize ...
sum.hist.flowering <- hist.flowering %>% group_by(year) %>% summarize(mdPrFL = median(prFL), lo50PrFL = quantile(prFL, 0.25), up50PrFL = quantile(prFL, 0.75), lo95PrFL = quantile(prFL, 0.025), up95PrFL = quantile(prFL, 0.975))

{cairo_pdf("output/figures/prFL-vs-time_linerange.pdf", width=6.5, height=3.5)

ggplot() + geom_linerange(data=sum.hist.flowering, aes(x=year, ymin=lo95PrFL, ymax=up95PrFL), color=park_palette("JoshuaTree")[5]) + geom_linerange(data=sum.hist.flowering, aes(x=year, ymin=lo50PrFL, ymax=up50PrFL), color=park_palette("JoshuaTree")[5], lwd=1) + geom_point(data=sum.hist.flowering, aes(x=year, y=mdPrFL), color="white", shape="-", size=2) +

geom_hline(yintercept=median(filter(hist.flowering, year<1930)$prFL), linetype=2) + geom_smooth(data=sub.hist.flowering, aes(x=year, y=prFL), method="loess", span=0.5, color=park_palette("JoshuaTree")[6], se=FALSE) + annotate("text", x=2021, y=1.05*median(filter(sub.hist.flowering, year<1930)$prFL), label="Median, 1900-1929", hjust=1, vjust=0, size=3) +

labs(x="Year", y="Probability of flowering") + dark_mode(theme_ipsum()) + theme(legend.position="none", plot.background=element_rect(color="black"))

}
dev.off()



