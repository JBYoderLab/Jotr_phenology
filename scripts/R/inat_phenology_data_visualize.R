# Scraping phenology-annotated iNat observations
# Assumes local environment 
# jby 2023.03.10

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

flow <- read.csv("output/flowering_obs_climate_v2_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)


# variant datasets -- dealing with the second flowering in 2019
flow2 <- flow %>% filter(!(year==2019.5 & flr==TRUE)) %>% mutate(year=floor(year)) # drop the late-flowering anomaly
flow3 <- flow
flow3$year[flow3$year==2019.5] <- 2019 # or merge 2019.5 into 2019?

glimpse(flow2) # 2,600 in our final working set

# split by subspecies
# swap input datasets to change --- current most trustworthy is flow2, ignoring 2019.5
yuja <- filter(flow2, type=="YUJA") 
yubr <- filter(flow2, type=="YUBR")

glimpse(yuja) # 1,124 obs (after iffy ones excluded)
glimpse(yubr) # 1,381 obs


#-------------------------------------------------------------------------
# visualize, if you like
inat_pheno_data <- read.csv("data/inat_phenology_data_subsp.csv", h=TRUE)

flr.raw.ln <- table(inat_pheno_data$year, inat_pheno_data$phenology) %>% as.data.frame() %>% rename(year=Var1, phenology=Var2, observations=Freq)

source("../shared/Rscripts/base_graphics.R")
library("ggdark")

if(!file.exists("output")) dir.create("output")
if(!file.exists("output/figures")) dir.create("output/figures")


# dark mode for presentations
{cairo_pdf("output/figures/Fig01E_iNat_obs_raw.pdf", width=5.25, height=3)

ggplot(flr.raw.ln, aes(x=year, y=observations, fill=phenology)) + geom_bar(stat="identity", position="dodge") + labs(x="Year of observation", y="iNaturalist records") + 

scale_fill_manual(values=park_palette("JoshuaTree")[c(6,1,3,5)], name="Phenology annotated") + 

theme_ipsum(base_size=10, axis_title_size=11, axis_title_just="rc") + theme(legend.position=c(0.2, 0.8), legend.key.size=unit(0.175,"in"), plot.margin=unit(c(0.1,0.1,0.05,0.1),"in"))

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
parks <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_parks_and_protected_lands_scale_rank")

rivers <- read_sf("../data/spatial/10m_physical/ne_10m_rivers_lake_centerlines_scale_rank", "ne_10m_rivers_lake_centerlines_scale_rank")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")
lakes <- read_sf("../data/spatial/10m_physical/ne_10m_lakes", "ne_10m_lakes")

sdm <- read_sf("../data/Yucca/jotr_BART_sdm_pres", "jotr_BART_sdm_pres")


statelabs <- data.frame(state=c("CA","NV","UT","AZ"), lon=c(-117.25,-116,-113.5,-113.45), lat=c(35.25,37.75,37.5,35.25))

# just the basics ....
{cairo_pdf(paste("output/figures/map_flr_obs.pdf", sep=""), width=6.5, height=6)

ggplot() + 

geom_sf(data=coast, color="slategray2", size=2.5) + 
geom_sf(data=states, fill="antiquewhite1", color=NA) + 
geom_text(data=statelabs, aes(label=state, x=lon, y=lat), color="white", size=12, alpha=0.75) +
geom_sf(data=urban, fill="antiquewhite2", color=NA) + 
geom_sf(data=states, fill=NA, color="antiquewhite4") + 
#geom_sf(data=parks, color=NA, fill=park_palette("JoshuaTree")[6], alpha=0.25) +

geom_sf(data=rivers, color="slategray3", size=0.4) +
geom_sf(data=lakes, fill="slategray3", color=NA) + 

geom_sf(data=sdm, fill=park_palette("JoshuaTree")[5], alpha=0.5, color=NA) +
#geom_sf(data=cole, fill=NA, color=park_palette("JoshuaTree")[5]) +

#geom_sf(data=parks, color=park_palette("JoshuaTree")[6], fill=NA, size=0.25) +

geom_point(data=flow2, aes(x=lon, y=lat, shape=flr, color=flr, size=flr), alpha=0.75) +

scale_shape_manual(values=c(21,20), labels=c("No evidence of flowering", "Buds, flowers, or fruit recorded"), name="Flowering activity") +
scale_size_manual(values=c(1.2,1), labels=c("No evidence of flowering", "Buds, flowers, or fruit recorded"), name="Flowering activity") +
scale_color_manual(values=park_palette("JoshuaTree")[c(5,3)], labels=c("No evidence of flowering", "Buds, flowers, or fruit recorded"), name="Flowering activity") + # nb these are false, true

annotation_scale(location = "br", width_hint = 0.3) + 
annotation_north_arrow(location = "br", which_north = "true", pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"), style = north_arrow_fancy_orienteering) +

coord_sf(xlim = c(-119, -112.75), ylim = c(33.25, 38), expand = FALSE) +

labs(x="Longitude", y="Latitude") +

theme_bw() + theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "slategray3"), axis.title=element_blank(), legend.position=c(0.81, 0.925), legend.box.margin=margin(0.01, 0.01, 0.01, 0.01, "in"), legend.key.height=unit(0.15, "in"))

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

