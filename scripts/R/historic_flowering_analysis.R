# Analyzing predicted historical flowering in Joshua tree
# Assumes local environment
# jby 2023.03.11

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
jotr.ri.files <- list.files("output/BART/RI.predictions", pattern=".bil", full=TRUE)

yubr.files <- list.files("output/BART/predictions.YUBR", pattern=".bil", full=TRUE)
yuja.files <- list.files("output/BART/predictions.YUJA", pattern=".bil", full=TRUE)

MojExt <- extent(-119, -112, 33, 38) # Mojave extent, maybe useful

# expert validation observations
vobs <- read.csv("data/Validation_obs_by_spp.csv")


# Jotr, without RI --------------------------------
jotr.histStack <- raster::stack(sapply(jotr.files, function(x) crop(raster::raster(x), MojExt)))
names(jotr.histStack) <- paste("prFL",1900:2022,sep=".")
projection(jotr.histStack)<-CRS("+init=epsg:4269")

jotr.histStack

jotr.maskHist <- mask(jotr.histStack, st_transform(sdm.pres[,2], crs=4269))

jotr.hist.flowering <- cbind(coordinates(jotr.maskHist), as.data.frame(jotr.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(jotr.hist.flowering)

# Jotr, with RI ------------------------------------
jotr.ri.histStack <- raster::stack(sapply(jotr.ri.files, function(x) crop(raster::raster(x), MojExt)))
names(jotr.ri.histStack) <- paste("prFL",1900:2022,sep=".")
projection(jotr.ri.histStack)<-CRS("+init=epsg:4269")

jotr.ri.histStack

jotr.ri.maskHist <- mask(jotr.ri.histStack, st_transform(sdm.pres[,2], crs=4269))

jotr.ri.hist.flowering <- cbind(coordinates(jotr.ri.maskHist), as.data.frame(jotr.ri.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="RI.prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(jotr.ri.hist.flowering)

jotr.flowering.histories <- jotr.hist.flowering |> left_join(jotr.ri.hist.flowering)

glimpse(jotr.flowering.histories)


write.table(jotr.flowering.histories, "output/historic_flowering_reconst_jotr.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

jotr.hist.flowering <- read.csv("output/historic_flowering_reconst_jotr.csv")

# comparing predictions with and without RI

table(jotr.flowering.histories$prFL>0.19, jotr.flowering.histories$RI.prFL>0.22) # hmmmm!
chisq.test(table(jotr.flowering.histories$prFL>0.19, jotr.flowering.histories$RI.prFL>0.22)) # p < 2.2e-16, natch

RIvNon <- jotr.flowering.histories |> group_by(year) |> do(broom::tidy(cor.test(~prFL+RI.prFL, data=., method="spearman"))) %>% ungroup()

glimpse(RIvNon)
RIvNon |> filter(p.value < 0.01, estimate < 0) # okay!
RIvNon |> filter(p.value < 0.01, estimate > 0) # whew! probabilities positively correlated in *all* years

ggplot(RIvNon, aes(x=estimate)) + geom_histogram() # and yeah the skew is in the direction we want


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
# jotr: 0.19
# jotr RI: 0.23

# YUBR: 0.30
# YUJA: 0.23

glimpse(vobs) 

valid <- NULL # initialize, in a lazy way

for(yr in unique(vobs$year)){

# yr <- 2005

vsub <- filter(vobs, year==yr)

vsub$jotr_prFlr <- raster::extract(jotr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$jotr_Flr <- raster::extract(jotr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.19

vsub$jotr_RI.prFlr <- raster::extract(jotr.ri.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$jotr_RI.Flr <- raster::extract(jotr.ri.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.23


#vsub$YUBR_prFlr <- raster::extract(yubr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
#vsub$YUBR_Flr <- raster::extract(yubr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.30

#vsub$YUJA_prFlr <- raster::extract(yuja.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
#vsub$YUJA_Flr <- raster::extract(yuja.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.23

valid <- rbind(vsub, valid)

}

glimpse(valid)
table(valid$obs_by)
head(valid) # cool
tail(valid)

write.table(valid, "output/expert_obs_validations.csv", sep=",", col.names=TRUE, row.names=FALSE)


# RANGE-WIDE -------------

# without RI
table(valid$obs_flowers, valid$jotr_Flr) # eesh
chisq.test(table(valid$obs_flowers, valid$jotr_Flr)) # n.s.

t.test(jotr_prFlr~obs_flowers, data=valid) # p = 0.08 hmmm

# jotr, with RI
table(valid$obs_flowers, valid$jotr_RI.Flr) # eesh
chisq.test(table(valid$obs_flowers, valid$jotr_RI.Flr)) # n.s.

t.test(jotr_RI.prFlr~obs_flowers, data=valid) # p = 0.1 hmmm

# BY SPECIES -------------
vYUBR <- filter(valid, type=="YUBR")
vYUJA <- filter(valid, type=="YUJA")

# without RI
t.test(jotr_prFlr~obs_flowers, data=vYUBR) # p = 8e-05 HEY
t.test(jotr_prFlr~obs_flowers, data=vYUJA) # n.s. and WRONG
chisq.test(table(vYUBR$obs_flowers, vYUBR$jotr_Flr)) # p = 0.02 (confirmed correct way)
chisq.test(table(vYUJA$obs_flowers, vYUJA$jotr_Flr)) # p = 0.0004 (WORSE than random wow)

# with RI
t.test(jotr_RI.prFlr~obs_flowers, data=vYUBR) # p = 0.001 OKAY
t.test(jotr_RI.prFlr~obs_flowers, data=vYUJA) # n.s. and WRONG
chisq.test(table(vYUBR$obs_flowers, vYUBR$jotr_RI.Flr)) # p = 0.09 (confirmed correct way)
chisq.test(table(vYUJA$obs_flowers, vYUJA$jotr_RI.Flr)) # n.s. whee


# can I visualize all this?

ggplot(valid, aes(y=jotr_prFlr, group=type, x=interaction(obs_flowers,type), shape=obs_flowers)) + geom_jitter() + geom_hline(yintercept=0.19, linetype=2) + scale_shape_manual(values=c(1,19)) 

ggplot(valid, aes(y=jotr_RI.prFlr, group=type, x=interaction(obs_flowers,type), shape=obs_flowers)) + geom_jitter() + geom_hline(yintercept=0.23, linetype=2) + scale_shape_manual(values=c(1,19)) 


#-------------------------------------------------------------------------
# observations versus same-year predictions

hist.flowering <- read.csv("output/historic_flowering_reconst_YUBR-YUJA.csv", h=TRUE) %>% mutate(type = factor(type, c("YUBR", "YUJA")))

glimpse(hist.flowering)

{cairo_pdf("output/figures/obs-vs-prediction_2021.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(jotr.hist.flowering, year==2021), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.5)) + 

geom_point(data=filter(obs, year==2021), aes(x=lon, y=lat, shape=flr, color=flr, size=flr)) + 
scale_shape_manual(values=c(21,20), name="Flowers seen?") + 
scale_size_manual(values=c(1.1,1), name="Flowers seen?") + 
scale_color_manual(values=park_palette("JoshuaTree")[c(5,3)], name="Flowers seen?") +

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()

{cairo_pdf("output/figures/obs-vs-prediction_2022.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(jotr.hist.flowering, year==2022), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.5)) + 

geom_point(data=filter(obs, year==2022), aes(x=lon, y=lat, shape=flr, color=flr, size=flr)) + 
scale_shape_manual(values=c(21,20), name="Flowers seen?") + 
scale_size_manual(values=c(1.1,1), name="Flowers seen?") + 
scale_color_manual(values=park_palette("JoshuaTree")[c(5,3)], name="Flowers seen?") +

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()





# all years with observations, global
{cairo_pdf("output/figures/obs-vs-prediction_Jotr.pdf", width=9.5, height=7.5)

ggplot() + 

geom_tile(data=filter(jotr.hist.flowering, year%in%2008:2022), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.5)) + 

geom_point(data=filter(obs, year%in%2008:2022), aes(x=lon, y=lat, shape=flr, color=flr), size=0.75) + 
scale_shape_manual(values=c(1,20), name="Flowers seen?") + 
scale_color_manual(values=c("gray60", park_palette("JoshuaTree")[2]), name="Flowers seen?") +

facet_wrap("year", nrow=3) + 

theme_minimal() + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.title=element_blank(), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6))

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
{cairo_pdf("output/figures/prediction_1982.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(jotr.hist.flowering, year==1982), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()



# year Joshua tree National Monument established
{cairo_pdf("output/figures/prediction_1936.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(jotr.hist.flowering, year==1936), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()

# earliest year in reconstruction
{cairo_pdf("output/figures/prediction_1900.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(jotr.hist.flowering, year==1900), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


#-------------------------------------------------------------------------
# TRENDS within cells

# cor(flr, yr) --- crude but it's a start -----------------

# not sure this is the best approach?
cellcors <- jotr.hist.flowering %>% group_by(lon,lat) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman"))) %>% ungroup() # rangewide model

glimpse(cellcors)

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
{cairo_pdf("output/figures/prFL-vs-time_map_jotr.pdf", width=6, height=5)

ggplot(cellcors, aes(x=lon, y=lat, fill=estimate)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="div", palette=1, direction=1, name=expression("Spearman's"~rho~", Pr(flowers) vs. time")) + labs(x="Longitude", y="Latitude") + 

theme_minimal() + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6))

}
dev.off()

# trends within single cells

declining <- cellcors %>% filter(p.value < 0.05, estimate < 0) 
increasing <- cellcors %>% filter(p.value < 0.01, estimate > 0) 
stable <- cellcors %>% filter(estimate > -0.05, estimate < 0.05) 


# map with example points highlighted
{cairo_pdf("output/figures/prFL-vs-time_map_jotr_examples.pdf", width=6, height=5)

ggplot(cellcors, aes(x=lon, y=lat, fill=estimate)) + geom_tile() + 

#coord_fixed() + 

geom_point(x=declining$lon[1], y=declining$lat[1], pch=21, color=park_palette("JoshuaTree")[2]) + 
geom_point(x=increasing$lon[1], y=increasing$lat[1], pch=21, color=park_palette("JoshuaTree")[2]) + 
geom_point(x=stable$lon[1], y=stable$lat[1], pch=21, color=park_palette("JoshuaTree")[2]) + 

scale_fill_distiller(type="div", palette=1, direction=1, name=expression("Spearman's"~rho~", Pr(flowers) vs. time")) + labs(x="Longitude", y="Latitude") + 

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


{cairo_pdf("output/figures/flowering_decline_example.pdf", width=4, height=3.5)

ggplot(filter(jotr.hist.flowering, lat == declining$lat[1], lon == declining$lon[1]), aes(x=year, y=prFL, color=prFL>0.24)) + geom_point() + geom_smooth(method="lm", color="#a6611a") +

scale_color_manual(values=c("white", "#80cdc1"), guide="none") + 

labs(x="Year", y="Predicted Pr(flowers)") +

dark_mode(theme_bw())

}
dev.off()

{cairo_pdf("output/figures/flowering_increase_example.pdf", width=4, height=3.5)

ggplot(filter(jotr.hist.flowering, lat == increasing$lat[1], lon == increasing$lon[1]), aes(x=year, y=prFL, color=prFL>0.24)) + geom_point() + geom_smooth(method="lm", color="#018571") +

scale_color_manual(values=c("white", "#80cdc1"), guide="none") + 

labs(x="Year", y="Predicted Pr(flowers)") +

dark_mode(theme_bw())

}
dev.off()

{cairo_pdf("output/figures/flowering_stable_example.pdf", width=4, height=3.5)

ggplot(filter(jotr.hist.flowering, lat == stable$lat[1], lon == stable$lon[1]), aes(x=year, y=prFL, color=prFL>0.24)) + geom_point() + geom_smooth(method="lm", color="white") +

scale_color_manual(values=c("white", "#80cdc1"), guide="none") + 

labs(x="Year", y="Predicted Pr(flowers)") +

dark_mode(theme_bw())

}
dev.off()

#-------------------------------------------------------------------------
# Flowering years

# best-power prediction thresholds from the original models
# jotr: 0.19
# jotr with RI: 0.23

# YUBR: 0.30
# YUJA: 0.23

flyrs.jotr <- jotr.hist.flowering %>% group_by(lon,lat) %>% 
summarize(
	flyrs_all=length(which(prFL>=0.22)), 
	flyrs_1900_1929=length(which(prFL>=0.22 & year<=1929)), 
	flyrs_1990_2019=length(which(prFL>=0.22 & year>=1990 & year<=2019)),
	flyrs_RI_all=length(which(RI.prFL>=0.23)), 
	flyrs_RI_1900_1929=length(which(RI.prFL>=0.23 & year<=1929)), 
	flyrs_RI_1990_2019=length(which(RI.prFL>=0.23 & year>=1990 & year<=2019)),
	) |> mutate(
	flyrs_change=flyrs_1990_2019-flyrs_1900_1929,
	flyrs_RI_change=flyrs_RI_1990_2019-flyrs_RI_1900_1929
	) # rangewide model
range(jotr.hist.flowering$year) # remember: 1900-2022
glimpse(flyrs.jotr)

write.table(flyrs.jotr, "output/jotr_reconstructed_flowering_years.csv", sep=",", col.names=TRUE, row.names=FALSE)


# without RI
quantile(flyrs.jotr$flyrs_all, c(0.025,0.5,0.975))
quantile(flyrs.jotr$flyrs_all, c(0.025,0.5,0.975))/123

quantile(flyrs.jotr$flyrs_1900_1929, c(0.025,0.5,0.975))/30 # hmmm
quantile(flyrs.jotr$flyrs_1990_2019, c(0.025,0.5,0.975))/30

# with RI
quantile(flyrs.jotr$flyrs_RI_all, c(0.025,0.5,0.975))
quantile(flyrs.jotr$flyrs_RI_all, c(0.025,0.5,0.975))/123

quantile(flyrs.jotr$flyrs_RI_1900_1929, c(0.025,0.5,0.975))/30 # hmmm
quantile(flyrs.jotr$flyrs_RI_1990_2019, c(0.025,0.5,0.975))/30


ggplot(flyrs.jotr, aes(x=flyrs_1900_1929)) + geom_histogram(bins=30)
ggplot(flyrs.jotr, aes(x=flyrs_1990_2019)) + geom_histogram(bins=30)


{cairo_pdf("output/figures/flowering-years_jotr.pdf", width=3.5, height=2.5)

ggplot(flyrs.jotr, aes(x=flyrs_all)) + geom_histogram(fill="#ccece6") + labs(x="Projected flowering years, 1900-2022", y="Grid cells") + 

geom_vline(xintercept=median(flyrs.jotr$flyrs_all), color="#006d2c") +

dark_mode(theme_minimal()) + theme(legend.position="none", plot.background=element_rect(color="black", fill="black"), plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()


# Flowering years, mapped
{cairo_pdf("output/figures/flyrs_map_jotr.pdf", width=6, height=5)

ggplot(flyrs.jotr, aes(x=lon, y=lat, fill=flyrs_all)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="seq", palette=2, direction=1, name="Years flowering predicted") + labs(x="Longitude", y="Latitude", title="Predicted flowering years, 1900-2022") + 

theme_minimal() + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6))

}
dev.off()

# Change in flowering years, mapped
{cairo_pdf("output/figures/flyrs_change_map_jotr.pdf", width=6, height=5)

ggplot(flyrs.jotr, aes(x=lon, y=lat, fill=flyrs_1990_2019-flyrs_1900_1929)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="div", palette=1, direction=1, name="Change in flowering years") + labs(x="Longitude", y="Latitude", title="Change in flowering years, 1990-2019 vs 1900-1929") + 

theme_minimal() + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6))

}
dev.off()

{cairo_pdf("output/figures/RI_flyrs_map_jotr.pdf", width=6, height=5)

ggplot(flyrs.jotr, aes(x=lon, y=lat, fill=flyrs_RI_all)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="seq", palette=2, direction=1, name="Years flowering predicted") + labs(x="Longitude", y="Latitude", title="Predicted flowering years, 1900-2022") + 

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

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


