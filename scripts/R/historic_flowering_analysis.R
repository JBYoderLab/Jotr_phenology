# Analyzing predicted historical flowering in Joshua tree
# Assumes local environment
# jby 2022.07.05

# starting up ------------------------------------------------------------

# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("rgdal")
library("raster")
library("sp")
library("rgeos")
library("sf")
library("ggspatial")

library("gdalUtilities")

library("rnaturalearth")
library("rnaturalearthdata")

library("mapproj")
library("maptools")
library("maps")
library("ggmap")

library("hexbin")

source("../shared/Rscripts/base.R") # my special mix of personal functions
source("../shared/Rscripts/base_graphics.R") # my special mix of personal functions

library("ggdark")

#-------------------------------------------------------------------------
# load up and organize historical data

# Jotr SDM
sdm.pres <- read_sf("../data/Yucca/jotr_BART_sdm_pres", "jotr_BART_sdm_pres")

# flowering observation data
obs <- read.csv("output/flowering_obs_rasterized.csv") %>% mutate(y2 = year) %>% mutate(year=floor(year))

# raster files of predicted prFL
files <- list.files("output/BART/predictions", pattern=".bil", full=TRUE)

MojExt <- extent(-119, -112, 33, 38)

histStack <- raster::stack(sapply(files, function(x) crop(raster::raster(x), MojExt)))
names(histStack) <- paste("prFL",1901:2022,sep=".")
projection(histStack)<-CRS("+init=epsg:4269")

histStack

maskHist <- mask(histStack, st_transform(sdm.pres, crs=4269))

hist.flowering <- cbind(coordinates(maskHist), as.data.frame(maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(hist.flowering)

write.table(hist.flowering, "output/historic_flowering_reconst.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

#-------------------------------------------------------------------------
# observations versus same-year predictions

{cairo_pdf("output/figures/obs-vs-prediction_2021.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==2021), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

geom_point(data=obs, aes(x=lon, y=lat, shape=flr, color=flr), size=0.5) + 
scale_shape_manual(values=c(1,20), name="Flowers observed") + 
scale_color_manual(values=c("gray60", park_palette("JoshuaTree")[2]), name="Flowers observed") +

coord_fixed() + labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


{cairo_pdf("output/figures/obs-vs-prediction.pdf", width=9.5, height=7)

ggplot() + 

geom_tile(data=filter(hist.flowering, year%in%2010:2021), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

geom_point(data=obs, aes(x=lon, y=lat, shape=flr, color=flr), size=0.5) + 
scale_shape_manual(values=c(1,20), name="Flowers observed") + 
scale_color_manual(values=c("gray60", park_palette("JoshuaTree")[2]), name="Flowers observed") +

coord_fixed() + facet_wrap("year") + 

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.title=element_blank(), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()

{cairo_pdf("output/figures/prediction_2022.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==2022), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

coord_fixed() + labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


{cairo_pdf("output/figures/prediction_1982.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==1982), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

coord_fixed() + labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


{cairo_pdf("output/figures/prediction_1910.pdf", width=6, height=5)

ggplot() + 

geom_tile(data=filter(hist.flowering, year==1910), aes(x=lon, y=lat, fill=prFL)) + scale_fill_distiller(type="seq", palette=2, direction=1, name="Pr(flowers)", limits=c(0,1), breaks=seq(0,1,0.2)) + 

coord_fixed() + labs(x="Latitude", y="Longitude") +

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

}
dev.off()


#-------------------------------------------------------------------------
# trends within cells

cellcors <- hist.flowering %>% group_by(lon,lat) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman")))

glimpse(cellcors)

{cairo_pdf("output/figures/prFL-vs-time_correlations.pdf", width=3.5, height=2.5)

ggplot(cellcors, aes(x=estimate, fill=p.value<0.01, color=p.value<0.01)) + geom_histogram(color=NA) + 

scale_fill_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +
scale_color_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +

annotate("text", x=0.1, y=150, label="All cells", color=park_palette("JoshuaTree")[7], fontface="bold", family="Arial Narrow") + 

annotate("text", x=-0.15, y=50, label="Cells with\ncorrelation\np < 0.01", color=park_palette("JoshuaTree")[2], fontface="bold", family="Arial Narrow", lineheight=0.8) + 

geom_vline(xintercept=0) +
labs(x=expression("Spearman's"~rho), y="Grid cells") + dark_mode(theme_minimal()) + theme(legend.position="none", plot.background=element_rect(color="black", fill="black"))

}
dev.off()

# Plot it on a map
{cairo_pdf("output/figures/prFL-vs-time_map.pdf", width=6, height=5)

ggplot(cellcors, aes(x=lon, y=lat, fill=estimate)) + geom_tile() + coord_fixed() + scale_fill_distiller(type="div", palette=1, direction=1, name=expression("Spearman's"~rho)) + labs(x="Longitude", y="Latitude", title="Correlation between probability of flowering and time") + 

dark_mode(theme_minimal()) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), plot.background=element_rect(color="black", fill="black"))

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

labs(x="Year", y="Probability of flowering") + theme_ipsum() + theme(legend.position="none")

}
dev.off()



