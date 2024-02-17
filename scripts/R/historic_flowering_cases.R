# Illustrating predicted flowering in Joshua tree for specific years
# Assumes local environment
# jby 2024.02.09

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("raster")
library("sp")
library("sf")
library("hexbin")
library("embarcadero")
library("rnaturalearth")
library("ggnewscale") # hmmm

source("../shared/Rscripts/base.R") # my special mix of personal functions
source("../shared/Rscripts/base_graphics.R") # my special mix of personal functions

library("ggdark")

#-------------------------------------------------------------------------
# load up and organize historical data

# Jotr predicted flowering

jotr.histStack <- raster::stack("output/BART/jotr_BART_predicted_flowering_1900-2023.grd")
jotr.ri.histStack <- raster::stack("output/BART/jotr_BART_RI_predicted_flowering_1900-2023.grd")

yubr.histStack <- raster::stack("output/BART/YUBR_BART_predicted_flowering_1900-2023.grd")
yuja.histStack <- raster::stack("output/BART/YUJA_BART_predicted_flowering_1900-2023.grd")

# check what we have here
plot(jotr.histStack[[1]])

# dataframe format, base and RI
jotr.hist.flowering <- read.csv("output/historic_flowering_reconst_jotr.csv")


# Jotr SDM and species boundaries
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range.shp")


#-------------------------------------------------------------------------
# Historical predictions

# Flowering in X year, mapped

# map elements
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range.shp")

parks <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_parks_and_protected_lands_scale_rank")

rivers <- read_sf("../data/spatial/10m_physical/ne_10m_rivers_lake_centerlines_scale_rank", "ne_10m_rivers_lake_centerlines_scale_rank")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")
lakes <- read_sf("../data/spatial/10m_physical/ne_10m_lakes", "ne_10m_lakes")


urban <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_urban_areas")
states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_admin_1_states_provinces")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")
JTNP <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_parks_and_protected_lands_scale_rank") %>% filter(unit_code=="JOTR")

topo <- projectRaster(raster("../data/spatial/10m_raster/US_MSR_10M/US_MSR.tif"), jotr.histStack[[1]])

topo_zoom <- crop(topo, extent(-117, -115, 33.5, 34.5))
topo_zoom

topo.df <- cbind(coordinates(raster(topo_zoom)), as.data.frame(raster(topo_zoom))) %>% rename(lon=x, lat=y, shade=layer) 


# Downey 1997 --- "Blooming Joshua trees greet visitors" --- Press-Enterprise
# Reynolds 2019 --- "Joshua Tree shows off its piece of the super bloom" --- LAT
# McKinney 1988 --- "Wildflowers, Joshua trees are bursting into bloom" --- LAT
# James 2013 --- "Record bloom at Joshua Tree" --- Desert Sun
# ?? 2005 --- "The Joshua Tree in Bloom, and Other Desert Flora." --- NYT


{cairo_pdf("output/figures/Historic_cases.pdf", width=4, height=7)

ggplot() + 
#geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
geom_sf(data=urban, fill="antiquewhite4", color=NA, lwd=0.1) + 

geom_tile(data=filter(jotr.hist.flowering, year%in%c(1987:1989, 1996:1998, 2004:2006, 2012:2014, 2018:2020)), aes(x=lon, y=lat, fill=prFL)) + 

facet_wrap("year", nrow=5) +

geom_sf(data=JTNP, fill=NA, color="white", linewidth=0.75) +

# water
geom_sf(data=rivers, color="slategray3", size=0.4) +
geom_sf(data=lakes[,-9], fill="slategray3", color=NA) + 

# labels
geom_text(data=data.frame(label="Joshua Tree NP", lon=-115.7, lat=34.2, year=1987), aes(label=label, x=lon, y=lat), color="white", size=3) +
geom_text(data=data.frame(label="Palm Springs", lon=-116, lat=33.6, year=1987), aes(label=label, x=lon, y=lat), color="antiquewhite4", size=3) +

scale_fill_distiller(type="seq", palette="Greens", direction=1, name="Pr(flowering)") + labs(x="Longitude", y="Latitude") + 

coord_sf(xlim = c(-116.5, -115.2), ylim = c(33.5, 34.3), expand = FALSE) +

theme_minimal(base_size=14) + theme(legend.position="bottom", legend.key.width=unit(0.35, "inches"), legend.key.height=unit(0.15, "in"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.01,0.01,0.01,0.01), "inches"), panel.spacing=unit(0.05,"inches"), legend.spacing.y=unit(0.05,"inches"), legend.box="horizontal", legend.box.spacing=unit(0.05, "inches"), legend.text=element_text(size=10), legend.title=element_text(size=12), strip.text=element_text(size=12), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())

}
dev.off()


# just the one 
ggplot() + 
#geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
geom_sf(data=urban, fill="antiquewhite4", color=NA, lwd=0.1) + 

geom_tile(data=filter(jotr.hist.flowering, year%in%c(1994:1996)), aes(x=lon, y=lat, fill=prFL)) + 

facet_wrap("year", nrow=5) +

geom_sf(data=JTNP, fill=NA, color="white", linewidth=0.75) +

# water
geom_sf(data=rivers, color="slategray3", size=0.4) +
geom_sf(data=lakes[,-9], fill="slategray3", color=NA) + 

# labels
geom_text(data=data.frame(label="Joshua Tree NP", lon=-115.7, lat=34.2, year=1987), aes(label=label, x=lon, y=lat), color="white", size=3) +
geom_text(data=data.frame(label="Palm Springs", lon=-116, lat=33.6, year=1987), aes(label=label, x=lon, y=lat), color="antiquewhite4", size=3) +

scale_fill_distiller(type="seq", palette="Greens", direction=1, name="Pr(flowering)") + labs(x="Longitude", y="Latitude") + 

coord_sf(xlim = c(-116.5, -115.2), ylim = c(33.5, 34.3), expand = FALSE) +

theme_minimal(base_size=14) + theme(legend.position="bottom", legend.key.width=unit(0.35, "inches"), legend.key.height=unit(0.15, "in"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.01,0.01,0.01,0.01), "inches"), panel.spacing=unit(0.05,"inches"), legend.spacing.y=unit(0.05,"inches"), legend.box="horizontal", legend.box.spacing=unit(0.05, "inches"), legend.text=element_text(size=10), legend.title=element_text(size=12), strip.text=element_text(size=12), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())

