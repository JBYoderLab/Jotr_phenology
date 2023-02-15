# Scraping phenology-annotated iNat observations
# Assumes local environment 
# jby 2023.01.04

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
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
non.y <- try(get_inat_obs(quality="research", place_id=53170, taxon_id=47785, term_id=12, term_value_id=21, year=y, maxresults=1e4))
Sys.sleep(5)  


if(class(bud.y)=="data.frame") bud.o <- bud.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Flower Budding", year=gsub("(\\d{4})-.+","\\1", observed_on)) else bud.o <- NULL

if(class(flo.y)=="data.frame") flo.o <- flo.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else flo.o <- NULL

if(class(fru.y)=="data.frame") fru.o <- fru.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="Fruiting", year=gsub("(\\d{4})-.+","\\1", observed_on)) else fru.o <- NULL

if(class(non.y)=="data.frame") non.o <- non.y %>% dplyr::select(scientific_name, latitude, longitude, url, image_url, observed_on) %>% mutate(phenology="No Evidence of Flowering", year=gsub("(\\d{4})-.+","\\1", observed_on)) else non.o <- NULL


inat_pheno_data <- rbind(inat_pheno_data, bud.o, flo.o, fru.o, non.o)

if(!file.exists("data")) dir.create("data")

write.table(inat_pheno_data, "data/inat_phenology_data.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)

} 

# expect error messages if searches return zero obs with a given phenology status; this may not be a problem, but see what the final data table looks like
glimpse(inat_pheno_data)
table(inat_pheno_data$year, inat_pheno_data$phenology)


#-------------------------------------------------------------------------
# separate (sub)species

# identify (sub)species
# inat_pheno_data <- read.csv("data/inat_phenology_data.csv", h=TRUE)
jtssps <- read_sf("data/Jotr_ssp_range.kml")

obs_in_ssp <- st_join(st_as_sf(inat_pheno_data, coords=c("longitude", "latitude"), crs=crs(jtssps)), jtssps, join = st_within) %>% cbind(inat_pheno_data[,c("longitude", "latitude")]) %>% as.data.frame(.) %>% dplyr::select(-geometry, -Description) %>% mutate(Name = gsub("(\\w+) Joshua tree range", "\\1", Name)) %>% rename(type=Name) %>% dplyr::select(latitude, longitude, type, year, observed_on, phenology)

glimpse(obs_in_ssp)
table(obs_in_ssp$phenology, obs_in_ssp$type)
table(obs_in_ssp$type)

write.table(obs_in_ssp, "data/inat_phenology_data_subsp.csv", sep=",", col.names=TRUE, row.names=FALSE)



# same treatment for CIS field observations
cis <- read.csv("data/CIS_obsv.csv") %>% mutate(lat = lat_deg+lat_min/60, lon = -(lon_deg+lon_min/60), year=as.numeric(gsub("\\d+/\\d+/(\\d+)", "\\1", date))) %>% dplyr::select(lat, lon, year, location, obs_flowers, obs_fruit, obs_moths, obs_no_flowers)

glimpse(cis)

cis_in_spp <- st_join(st_as_sf(cis, coords=c("lon", "lat"), crs=crs(jtssps)), jtssps, join = st_within) %>% cbind(cis[,c("lon", "lat")]) %>% as.data.frame(.) %>% dplyr::select(-geometry, -Description) %>% mutate(Name = gsub("(\\w+) Joshua tree range", "\\1", Name)) %>% rename(type=Name) %>% dplyr::select(lat, lon, location, type, year, obs_flowers, obs_fruit, obs_moths, obs_no_flowers)

glimpse(cis_in_spp)

write.table(cis_in_spp, "data/CIS_obs_by_spp.csv", col.names=TRUE, row.names=FALSE, sep=",")

# rasterize!

valid <- read.csv("data/Validation_obs_by_spp.csv")
glimpse(valid)

valid.rst <- data.frame(matrix(0,0,4))
names(valid.rst) <- c("lon","lat","year","flr")

prism_temp_rast <- raster("data/PRISM/annual/ppt_Mojave_2010Q1.bil")

# then ...
for(yr in 1:length(unique(valid$year))){

if(nrow(dplyr::filter(valid, year==unique(valid$year)[yr], ((obs_by=="CIS" & (obs_flowers|obs_fruit)) | (obs_by=="Ray Yeager" & obs_flowers)) ))> 0){
yes <- 2*(rasterize(dplyr::filter(valid, year==unique(valid$year)[yr], ((obs_by=="CIS" & (obs_flowers|obs_fruit)) | (obs_by=="Ray Yeager" & obs_flowers)) )[,c("lon","lat")], prism_temp_rast, fun=sum, background=0) > 0)

yearyes <- rasterToPoints(yes+no, fun=function(x){x>1})
}else{
yearyes <- NULL
}


if(nrow(dplyr::filter(valid, year==unique(valid$year)[yr], !((obs_by=="CIS" & (obs_flowers|obs_fruit)) | (obs_by=="Ray Yeager" & obs_flowers)) ))> 0){
no <- rasterize(dplyr::filter(valid, year==unique(valid$year)[yr], !((obs_by=="CIS" & (obs_flowers|obs_fruit)) | (obs_by=="Ray Yeager" & obs_flowers)) )[,c("lon","lat")], prism_temp_rast, fun=sum, background=0) > 0

yearno <- rasterToPoints(yes+no, fun=function(x){x==1})

}else{
yearno <- NULL
}

valid.rst <- rbind(valid.rst, data.frame(lon=c(yearyes[,"x"],yearno[,"x"]), lat=c(yearyes[,"y"],yearno[,"y"]), year=unique(valid$year)[yr], flr=rep(c(TRUE,FALSE),c(nrow(yearyes),nrow(yearno)))))# need to solve this yet

}

valid.rst.ssp <- st_join(st_as_sf(valid.rst, coords=c("lon", "lat"), crs=crs(jtssps)), jtssps, join = st_within) %>% cbind(valid.rst[,c("lon", "lat")]) %>% as.data.frame(.) %>% dplyr::select(-geometry, -Description) %>% rename(type=Name) %>% dplyr::select(lon, lat, type, year, flr)

glimpse(valid.rst.ssp) # okay okay okay!
table(valid.rst.ssp$type, useNA="ifany")
table(valid.rst.ssp$year, valid.rst.ssp$flr) # think that looks good ...


glimpse(valid.rst) # okay, let's see how we feel about this




#-------------------------------------------------------------------------
# visualize, if you like
inat_pheno_data <- read.csv("data/inat_phenology_data.csv", h=TRUE)

flr.raw.ln <- table(inat_pheno_data$year, inat_pheno_data$phenology) %>% as.data.frame() %>% rename(year=Var1, phenology=Var2, observations=Freq)

source("../shared/Rscripts/base_graphics.R")
library("ggdark")

if(!file.exists("output")) dir.create("output")
if(!file.exists("output/figures")) dir.create("output/figures")


{cairo_pdf("output/figures/iNat_obs_raw.pdf", width=11, height=5)

ggplot(flr.raw.ln, aes(x=year, y=observations, fill=phenology)) + geom_bar(stat="identity", position="dodge") + labs(x="Year of observation", y="iNat records (research grade)") + scale_fill_manual(values=park_palette("JoshuaTree")[c(6,7,3,5)], name="Phenology") + dark_mode(theme_ipsum(base_size=14, axis_title_size=20)) + theme(plot.background=element_rect(color="black"))

}
dev.off()

