# Scraping phenology-annotated iNat observations
# Assumes local environment 
# jby 2022.08.08

# starting up ------------------------------------------------------------

# setwd("/Volumes/GoogleDrive/Other computers/My MacBook Pro 2020/Documents/Academic/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("sf")
library("raster")

source("scripts/R/get_inat.R") # attaches rinat and hacks the key function

#-------------------------------------------------------------------------
# Pull down iNat observations of Joshua tree with specific phenology code

# term id: 12 for Plant Phenology then term_id_value: 13 =Flowering, 14 =Fruiting, 15 =Flower Budding

# trial run, to make sure it works as expected
test <- get_inat_obs(quality="research", place_id=53170, taxon_id=47785, term_id=12, term_value_id=13, year=2021, maxresults=1e4) 

glimpse(test)


# ACTUALLY RUN THE THING
inat_pheno_data <- data.frame(matrix(0,0,7))
names(inat_pheno_data) <- c("scientific_name", "latitude", "longitude", "url", "image_url", "observed_on", "phenology")

years <- 2010:2022 # to run everything

# to read back in and continue
# inat_pheno_data <- read.csv("data/inat_phenology_data.csv", h=TRUE)

# okay let's pull this stuff down already
# n.b. for-looping this borks up in a way that makes me suspect it's overloading the API
for(y in years){

# y <- 2022

bud.y <- try(get_inat_obs(quality="research", place_id=53170, taxon_id=47785, term_id=12, term_value_id=15, year=y, maxresults=1e4))
Sys.sleep(5) # throttling under the API limit, maybe?
flo.y <- try(get_inat_obs(quality="research", place_id=53170, taxon_id=47785, term_id=12, term_value_id=13, year=y, maxresults=1e4))
Sys.sleep(5)  
fru.y <- try(get_inat_obs(quality="research", place_id=53170, taxon_id=47785, term_id=12, term_value_id=14, year=y, maxresults=1e4))
Sys.sleep(5)  
non.y <- try(get_inat_obs(quality="research", place_id=53170, taxon_id=47785, term_id=12, without_term_value_id="13,14,15", year=y, maxresults=1e4))
Sys.sleep(5)  


if(class(bud.y)=="data.frame") bud.o <- bud.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Flower Budding", year=gsub("(\\d{4})-.+","\\1", observed_on)) else bud.o <- NULL

if(class(flo.y)=="data.frame") flo.o <- flo.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else flo.o <- NULL

if(class(fru.y)=="data.frame") fru.o <- fru.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Fruiting", year=gsub("(\\d{4})-.+","\\1", observed_on)) else fru.o <- NULL

if(class(non.y)=="data.frame") non.o <- non.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="No Evidence of Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else non.o <- NULL


inat_pheno_data <- rbind(inat_pheno_data, bud.o, flo.o, fru.o, non.o)

write.table(inat_pheno_data, "data/inat_phenology_data.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

} 

# expect error messages if searches return zero obs with a given phenology status; this may not be a problem, but see what the final data table looks like
glimpse(inat_pheno_data)
table(inat_pheno_data$year, inat_pheno_data$phenology)


#-------------------------------------------------------------------------
# separate (sub)species

# identify (sub)species
# flr.clim <- read.csv("output/flowering_obs_climate_v2.csv", h=TRUE)
jtssps <- read_sf("data/Jotr_range.kml")

obs_in_ssp <- st_join(st_as_sf(inat_pheno_data, coords=c("longitude", "latitude"), crs=crs(jtssps)), jtssps, join = st_within) %>% cbind(inat_pheno_data[,c("longitude", "latitude")]) %>% as.data.frame(.) %>% dplyr::select(-geometry, -Description) %>% mutate(Name = gsub("(\\w+) Joshua tree range", "\\1", Name)) %>% rename(type=Name) %>% dplyr::select(latitude, longitude, type, year, observed_on, phenology)

glimpse(obs_in_ssp)

write.table(obs_in_ssp, "data/inat_phenology_data_subsp.csv", sep=",", col.names=TRUE, row.names=FALSE)

table(obs_in_ssp$phenology, obs_in_ssp$type)
table(obs_in_ssp$type)


#-------------------------------------------------------------------------
# visualize, if you like
inat_pheno_data <- read.csv("data/inat_phenology_data.csv", h=TRUE)

flr.raw.ln <- table(inat_pheno_data$year, inat_pheno_data$phenology) %>% as.data.frame() %>% rename(year=Var1, phenology=Var2, observations=Freq)

source("../shared/Rscripts/base_graphics.R")
library("ggdark")

{cairo_pdf("output/figures/iNat_obs_raw.pdf", width=11, height=5)

ggplot(flr.raw.ln, aes(x=year, y=observations, fill=phenology)) + geom_bar(stat="identity", position="dodge") + labs(x="Year of observation", y="iNat records (research grade)") + scale_fill_manual(values=park_palette("JoshuaTree")[c(6,7,3,5)], name="Phenology") + dark_mode(theme_ipsum(base_size=14, axis_title_size=20)) + theme(plot.background=element_rect(color="black"))

}
dev.off()

