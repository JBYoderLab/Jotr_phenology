# Scraping phenology-annotated iNat observations
# Assumes local environment 
# jby 2024.01.24

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("sf")
library("raster")
library("ggspatial")

source("../shared/Rscripts/base_graphics.R")

#-------------------------------------------------------------------------
# load data

flow <- read.csv("output/flowering_obs_climate_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)


# variant datasets -- dealing with the second flowering in 2018
flow2 <- flow %>% filter(!(year==2018.5 & flr==TRUE), year>=2008) %>% mutate(year=floor(year)) # drop the late-flowering anomaly
flow3 <- flow |> filter(year>=2008)
flow3$year[flow3$year==2018.5] <- 2018 # or merge 2018.5 into 2018?

glimpse(flow2) # 2,632 in our final working set (more than this with 2023!)

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="YUJA") 
yubr <- filter(flow2, type=="YUBR")

glimpse(yuja) # 1,172 obs (after iffy ones excluded)
glimpse(yubr) # 1,460 obs


#-------------------------------------------------------------------------
# visualize all phenology annotations
inat_pheno_data <- read.csv("data/inat_phenology_data_subsp.csv", h=TRUE) |> filter(year>=2008)

table(inat_pheno_data$year)

flr.raw.ln <- table(inat_pheno_data$year, inat_pheno_data$phenology) %>% as.data.frame() %>% rename(year=Var1, phenology=Var2, observations=Freq) %>% mutate(year=as.numeric(as.character(year)))
glimpse(flr.raw.ln)

source("../shared/Rscripts/base_graphics.R")
library("ggdark")

if(!file.exists("output")) dir.create("output")
if(!file.exists("output/figures")) dir.create("output/figures")


# dark mode for presentations
{cairo_pdf("output/figures/Fig01E_iNat_obs_raw.pdf", width=5.25, height=3)

ggplot(flr.raw.ln, aes(x=year, y=observations, fill=phenology)) + geom_bar(stat="identity", position="dodge") + labs(x="Year of observation", y="iNaturalist records") +

scale_fill_manual(values=park_palette("JoshuaTree")[c(6,1,3,5)], name="Phenology annotated") + 

theme_minimal(base_size=13) + theme(legend.position=c(0.25, 0.8), legend.key.size=unit(0.15,"in"), plot.margin=unit(c(0.1,0.1,0.05,0.1),"in"), axis.text.x = element_text(angle=45))

}
dev.off()

{png("output/figures/iNat_obs_example.png", width=565, height=540)

ggplot(filter(flr.raw.ln, year>2014), aes(x=factor(year), y=observations, fill=phenology)) + geom_bar(stat="identity", position="dodge") + labs(x="Year of observation", y="iNaturalist records") +

scale_fill_manual(values=park_palette("JoshuaTree")[c(6,1,3,5)], name="Phenology annotated") + 

theme_minimal(base_size=26) + theme(legend.position=c(0.35, 0.8), legend.key.size=unit(0.2,"in"), plot.margin=unit(c(0.1,0.1,0.05,0.25),"in"), axis.text.x = element_text(angle=45))

}
dev.off()

#-------------------------------------------------------------------------
# visualize predictor differences
# plot the predictors in raw observations ...

# Full range
jotr.preds <- c("tmaxW0vW1", "pptY0", "tmaxW0", "tminW0") # output from varimp(), modeling script

jotr.pred.plot <- flow2 %>% dplyr::select(year, flr, all_of(jotr.preds)) %>% pivot_longer(all_of(jotr.preds), names_to="Predictor", values_to="Value")


{cairo_pdf(file="output/figures/RImod_best_predictors_Jotr.pdf", width=6, height=3.5)

ggplot(jotr.pred.plot, aes(x=flr, y=Value, color=Predictor)) + geom_jitter(alpha=0.1, size=0.25) + geom_boxplot(alpha=0.5, aes(fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,2,2,2)], guide="none") +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,2,2,2)], guide="none") +

scale_x_discrete(labels=c("N","Y")) + 

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + theme_bw(base_size=15) + theme()


}
dev.off()

# dark mode for presentations
{cairo_pdf(file="output/figures/RImod_best_predictors_Jotr_present.pdf", width=6, height=3.5)

ggplot(jotr.pred.plot, aes(x=flr, y=Value, color=Predictor)) + geom_jitter(alpha=0.1, size=0.25) + geom_boxplot(alpha=0.5, aes(fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,2,2,2)], guide="none") +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,2,2,2)], guide="none") +

scale_x_discrete(labels=c("N","Y")) + 

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + dark_mode(theme_ipsum(base_size=12, axis_title_size=14, strip_text_size=10)) + theme(plot.background=element_rect(color="black"))


}
dev.off()



# YUBR only
yubr.preds <- c("tmaxW0vW1", "pptY0", "tmaxW0", "tminW0") # output from varimp(), modeling script

yubr.pred.plot <- yubr %>% dplyr::select(year, flr, all_of(yubr.preds)) %>% pivot_longer(all_of(yubr.preds), names_to="Predictor", values_to="Value")

# dark mode for presentations
{cairo_pdf(file="output/figures/RImod_best_predictors_yubr_present.pdf", width=6, height=3.5)

ggplot(yubr.pred.plot, aes(x=flr, y=Value, color=Predictor)) + geom_jitter(alpha=0.1, size=0.25) + geom_boxplot(alpha=0.5, aes(fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,2,2,2)], guide="none") +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,2,2,2)], guide="none") +

scale_x_discrete(labels=c("N","Y")) + 

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + dark_mode(theme_ipsum(base_size=12, axis_title_size=14, strip_text_size=10)) + theme(plot.background=element_rect(color="black"))


}
dev.off()


# YUJA only
yuja.preds <- c("vpdmaxW0vW1", "tmaxW0", "vpdmaxW0", "pptW0", "pptW0W1") # output from varimp(), modeling script

yuja.pred.plot <- flow2 %>% dplyr::select(year, flr, all_of(yuja.preds)) %>% pivot_longer(all_of(yuja.preds), names_to="Predictor", values_to="Value")

# dark mode for presentations
{cairo_pdf(file="output/figures/RImod_best_predictors_yuja_present.pdf", width=7, height=3.5)

ggplot(yuja.pred.plot, aes(x=flr, y=Value, color=Predictor)) + geom_jitter(alpha=0.1, size=0.25) + geom_boxplot(alpha=0.5, aes(fill=Predictor), width=0.5) + 

scale_color_manual(values=park_palette("JoshuaTree")[c(1,1,2,6,6)], guide="none") +
scale_fill_manual(values=park_palette("JoshuaTree")[c(1,1,2,6,6)], guide="none") +

scale_x_discrete(labels=c("N","Y")) + 

facet_wrap("Predictor", nrow=1, scale="free_y") + labs(x="Flowers observed?", y="Predictor value") + dark_mode(theme_ipsum(base_size=10, axis_title_size=14, strip_text_size=10)) + theme(plot.background=element_rect(color="black"), panel.spacing.x=unit(0.15, "inches"))


}
dev.off()


#-------------------------------------------------------------------------
# Map observations


# geography stuff ...................
states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_admin_1_states_provinces")
urban <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_urban_areas")
parks <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_parks_and_protected_lands_scale_rank") |> filter(unit_name=="Joshua Tree NP")

rivers <- read_sf("../data/spatial/10m_physical/ne_10m_rivers_lake_centerlines_scale_rank", "ne_10m_rivers_lake_centerlines_scale_rank")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")
lakes <- read_sf("../data/spatial/10m_physical/ne_10m_lakes", "ne_10m_lakes")

sdm <- read_sf("../data/Yucca/Jotr_SDM2023_range.shp")


statelabs <- data.frame(state=c("CA","NV","UT","AZ"), lon=c(-116.75,-115.75,-113.5,-113.45), lat=c(35.25,37.65,37.5,35.5))
citylabs <- data.frame(name=c("Los Angeles"), lon=c(-118), lat=c(34))
parklabs <- data.frame(name=c("Joshua Tree NP"), lon=c(-115.1), lat=c(34.2))


# just the basics ....
{cairo_pdf(paste("output/figures/map_flr_obs.pdf", sep=""), width=6.5, height=6)

ggplot() + 

geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
geom_sf(data=states, fill="antiquewhite1", color=NA) + 
geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="antiquewhite3", size=16, alpha=0.75) +
geom_sf(data=urban, fill="antiquewhite2", color=NA) + 
geom_sf(data=states, fill=NA, color="antiquewhite4") + 
geom_sf(data=parks, color=park_palette("JoshuaTree")[5], fill=NA, size=1.5) +

geom_sf(data=rivers, color="slategray3", linewidth=0.4) +
geom_sf(data=lakes, fill="slategray3", color=NA) + 

geom_text(data=citylabs, aes(label=name, x=lon, y=lat), color="black", size=5.5, alpha=0.75) +
geom_text(data=parklabs, aes(label=name, x=lon, y=lat), color=park_palette("JoshuaTree")[5], size=5, alpha=0.75) +


geom_sf(data=sdm, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
#geom_sf(data=cole, fill=NA, color=park_palette("JoshuaTree")[5]) +

geom_point(data=flow2, aes(x=lon, y=lat, shape=flr, color=flr, size=flr), alpha=0.75) +

scale_shape_manual(values=c(21,20), labels=c("No flowering", "Buds, flowers, or fruit"), name="Flowering, 2008-22") +
scale_size_manual(values=c(1.3,1.1), labels=c("No flowering", "Buds, flowers, or fruit"), name="Flowering, 2008-22") +
scale_color_manual(values=c("#ffffcc", "#253494"), labels=c("No flowering", "Buds, flowers, or fruit"), name="Flowering, 2008-22") + # nb these are false, true

annotation_scale(location = "br", width_hint = 0.25, height = unit(0.15, "in"), text_cex=1) + 
annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering(text_size=16), height = unit(1,"in"), width = unit(0.75, "in")) +

coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38.25), expand = TRUE) +

labs(x="Longitude", y="Latitude") +

theme_bw(base_size=14) + theme(panel.grid=element_blank(), panel.background = element_rect(fill = "slategray3", linewidth=0.25), axis.title = element_blank(), axis.text=element_blank(), axis.ticks = element_blank(), legend.position=c(0.825, 0.92), legend.box.margin = margin(0.01, 0.01, 0.01, 0.01, "in"), legend.key.size=unit(0.125, "in"), legend.key = element_rect(fill="#B0AC85"), legend.background = element_rect(color="black", linewidth=0.1), plot.margin = margin(10,10,10,10))

}
dev.off()

# (sub)species ranges
spp.ranges <- read_sf("data/Jotr_ssp_range.kml")

yubr.pres <- st_zm(st_intersection(sdm.pres, st_transform(filter(spp.ranges, Name=="YUBR"), crs=crs(sdm.pres))))
yuja.pres <- st_zm(st_intersection(sdm.pres, st_transform(filter(spp.ranges, Name=="YUJA"), crs=crs(sdm.pres))))

{cairo_pdf(paste("output/figures/map_flr_obs_ssp.pdf", sep=""), width=6.5, height=6)

ggplot() + 

geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
geom_sf(data=states, fill="antiquewhite1", color=NA) + 
geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="antiquewhite3", size=16, alpha=0.75) +
geom_sf(data=urban, fill="antiquewhite2", color=NA) + 
geom_sf(data=states, fill=NA, color="antiquewhite4") + 
geom_sf(data=parks, color=park_palette("JoshuaTree")[5], fill=NA, size=1.5) +

geom_sf(data=rivers, color="slategray3", linewidth=0.4) +
geom_sf(data=lakes, fill="slategray3", color=NA) + 

geom_text(data=citylabs, aes(label=name, x=lon, y=lat), color="black", size=5.5, alpha=0.75) +
geom_text(data=parklabs, aes(label=name, x=lon, y=lat), color=park_palette("JoshuaTree")[5], size=5, alpha=0.75) +

geom_sf(data=yubr.pres, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
geom_sf(data=yuja.pres, fill=park_palette("JoshuaTree")[6], alpha=0.5, color=NA) +

#geom_sf(data=cole, fill=NA, color=park_palette("JoshuaTree")[5]) +

geom_point(data=flow2, aes(x=lon, y=lat, shape=flr, color=flr, size=flr), alpha=0.75) +

scale_shape_manual(values=c(21,20), labels=c("No flowering", "Buds, flowers, or fruit"), name="Flowering, 2008-22") +
scale_size_manual(values=c(1.3,1.1), labels=c("No flowering", "Buds, flowers, or fruit"), name="Flowering, 2008-22") +
scale_color_manual(values=c("#ffffcc", "#253494"), labels=c("No flowering", "Buds, flowers, or fruit"), name="Flowering, 2008-22") + # nb these are false, true

annotation_scale(location = "br", width_hint = 0.25, height = unit(0.15, "in"), text_cex=1) + 
annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering(text_size=16), height = unit(1,"in"), width = unit(0.75, "in")) +

coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38), expand = FALSE) +

labs(x="Longitude", y="Latitude") +

theme_bw(base_size=14) + theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "slategray3", linewidth=0.25), axis.title = element_blank(), axis.ticks = element_line(linewidth=0.2), legend.position=c(0.83, 0.925), legend.box.margin = margin(0.01, 0.01, 0.01, 0.01, "in"), legend.key.size=unit(0.15, "in"), legend.key = element_rect(fill="#B0AC85"), legend.background = element_rect(color="black", linewidth=0.1), plot.margin = margin(10,10,25,10))

}
dev.off()

# background map
{cairo_pdf(paste("output/figures/map_jotr_range.pdf", sep=""), width=6.5, height=6)

ggplot() + 

geom_sf(data=coast, color="slategray2", size=2.5) + 
geom_sf(data=states, fill="antiquewhite1", color=NA) + 
geom_text(data=statelabs, aes(label=state, x=lon, y=lat), color="white", size=12, alpha=0.75) +
geom_sf(data=urban, fill="antiquewhite2", color=NA) + 
geom_sf(data=states, fill=NA, color="antiquewhite4") + 
geom_sf(data=parks, color=NA, fill=park_palette("JoshuaTree")[6], alpha=0.25) +

geom_sf(data=rivers, color="slategray3", size=0.4) +
geom_sf(data=lakes, fill="slategray3", color=NA) + 

geom_sf(data=sdm, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
#geom_sf(data=cole, fill=NA, color=park_palette("JoshuaTree")[5]) +

geom_sf(data=parks, color=park_palette("JoshuaTree")[6], fill=NA, size=0.25) +


annotation_scale(location = "br", width_hint = 0.3) + 
annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering) +

coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38), expand = FALSE) +

labs(x="Longitude", y="Latitude") +

dark_mode(theme_bw()) + theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "slategray3"), plot.background = element_rect(color="black"))

}
dev.off()

#-------------------------------------------------------------------------
# Figures for a fruiting model

inat_pheno_data <- read.csv("data/inat_phenology_data_subsp.csv", h=TRUE) |> filter(year>=2008)
glimpse(inat_pheno_data)

table(inat_pheno_data$year)

# first, all phenology states
flr.raw.ln <- table(inat_pheno_data$year, inat_pheno_data$phenology) %>% as.data.frame() %>% rename(year=Var1, phenology=Var2, observations=Freq) %>% mutate(year=as.numeric(as.character(year)))
glimpse(flr.raw.ln)


allpheno <- ggplot(flr.raw.ln, aes(x=year, y=observations, fill=phenology)) + 
			geom_bar(stat="identity", position="stack") + labs(x="Year of observation", y="iNaturalist records") +

			scale_fill_manual(values=park_palette("JoshuaTree")[c(2,1,3,5)], name="Phenology", labels=c("Buds", "Flowers", "Fruit", "No flowering")) + 

			theme_minimal(base_size=14) + theme(legend.position=c(0.3, 0.75), legend.key.size=unit(0.2,"in"), plot.margin=unit(c(0.1,0.01,0.05,0.1),"in"), axis.text.x = element_blank(), axis.title.x = element_blank())

allpheno 


# then, proportion fruiting/flowering only
fruit.ln <- inat_pheno_data %>% filter(phenology!="No Evidence of Flowering") %>% group_by(year, phenology) %>% 
			summarize(N = n()) %>% pivot_wider(names_from=phenology, values_from=N) %>% 
			rename(Budding=`Flower Budding`) %>% replace_na(list(Fruiting=0, Budding=0)) %>% 
			mutate(N = Flowering + Fruiting + Budding, prop_fruit=Fruiting/N, prop_not_fruit=(Flowering+Budding)/N) %>% 
			dplyr::select(year, prop_fruit, prop_not_fruit) %>% pivot_longer(starts_with("prop"), names_to="phenology", values_to="prop")
fruit.ln # there we go, oy

fruitplot <- ggplot(fruit.ln, aes(x=year, y=prop, fill=phenology)) +
			 geom_bar(stat="identity", position="stack") + labs(x="Year of observation", y="Proportion of records") +
			 
			 scale_fill_manual(values=park_palette("JoshuaTree")[c(1,3)], guide="none") +
			 
			 annotate("text", label=c("Fruit", "Buds or flowers"), x=2008, y=c(0.01,0.51), color="white", angle=90, hjust=0, fontface="bold", size=3) +
			 
			 theme_minimal(base_size=14) + theme(legend.position=c(0.35, 0.8), legend.key.size=unit(0.2,"in"), plot.margin=unit(c(0.1,0.01,0.1,0.1),"in"), axis.text.x = element_text(angle=45))

fruitplot


# and finally the map

# geography stuff ...................
states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_admin_1_states_provinces")
urban <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_urban_areas")
parks <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_parks_and_protected_lands_scale_rank") |> filter(unit_name=="Joshua Tree NP")

rivers <- read_sf("../data/spatial/10m_physical/ne_10m_rivers_lake_centerlines_scale_rank", "ne_10m_rivers_lake_centerlines_scale_rank")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")
lakes <- read_sf("../data/spatial/10m_physical/ne_10m_lakes", "ne_10m_lakes")

sdm <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range.shp")


statelabs <- data.frame(state=c("CA","NV","UT","AZ"), lon=c(-116.75,-115.75,-113.5,-113.45), lat=c(35.25,37.65,37.5,35.5))
citylabs <- data.frame(name=c("Los Angeles"), lon=c(-118), lat=c(34))
parklabs <- data.frame(name=c("Joshua Tree NP"), lon=c(-115.1), lat=c(34.2))


# just the basics ....
fruitmap <- ggplot() + 

			geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
			geom_sf(data=states, fill="antiquewhite1", color=NA) + 
			geom_text(data=filter(statelabs, state!="UT"), aes(label=state, x=lon, y=lat), color="antiquewhite3", size=16, alpha=0.75) +
			geom_sf(data=states, fill=NA, color="antiquewhite4") + 

			geom_sf(data=rivers, color="slategray3", linewidth=0.4) +
			geom_sf(data=lakes, fill="slategray3", color=NA) + 

			geom_sf(data=sdm, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
			
			geom_point(data=filter(inat_pheno_data, phenology!="No Evidence of Flowering"), aes(x=longitude, y=latitude, shape=phenology=="Fruiting", color=phenology=="Fruiting", size=phenology=="Fruiting"), alpha=0.75) +

			scale_shape_manual(values=c(21,20), labels=c("Buds or flowers", "Fruit"), name="Fruiting, 2008-23") +
			scale_size_manual(values=c(1.7,1.5), labels=c("Buds or flowers", "Fruit"), name="Fruiting, 2008-23") +
			scale_color_manual(values=park_palette("JoshuaTree")[c(1,3)], labels=c("Buds or flowers", "Fruit"), name="Fruiting, 2008-23") + # nb these are false, true

			annotation_scale(location = "br", width_hint = 0.25, height = unit(0.15, "in"), text_cex=1) + 
			annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering(text_size=16), height = unit(1,"in"), width = unit(0.75, "in")) +

			coord_sf(xlim = c(-118.75, -113), ylim = c(33.25, 38.25), expand = TRUE) +

			labs(x="Longitude", y="Latitude") +

			theme_bw(base_size=14) + theme(panel.grid=element_blank(), panel.background = element_rect(fill = "slategray3", linewidth=0.25), axis.title = element_blank(), axis.text=element_blank(), axis.ticks = element_blank(), legend.position=c(0.825, 0.92), legend.box.margin = margin(0.01, 0.01, 0.01, 0.01, "in"), legend.key.size=unit(0.125, "in"), legend.key = element_rect(fill="#B0AC85"), legend.background = element_rect(color="black", linewidth=0.1), plot.margin = margin(10,10,10,10))

fruitmap


# big composite figure			
{cairo_pdf(paste("output/figures/fruiting_pheno_data.pdf", sep=""), width=11, height=6)

ggdraw() + draw_plot(allpheno, 0, 0.5, 0.47, 0.5) + draw_plot(fruitplot, 0, 0, 0.47, 0.5) + draw_plot(fruitmap, 0.45, 0, 0.55, 1) + draw_plot_label(label=c("A", "B", "C"), x=c(0,0,0.49), y=c(0.975,0.52,0.975), size=20)

}
dev.off()
