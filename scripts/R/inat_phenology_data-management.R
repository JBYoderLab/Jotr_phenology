# working with phenology-annotated iNat observations
# Assumes MAJEL environment 
# jby 2022.11.08

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("lubridate")

library("raster")
library("sf")

# Mojave crop extent (deliberately generous)
MojExt <- extent(-119, -112, 33, 38)

#-------------------------------------------------------------------------
# read in iNat observations compiled using inat_phenology_obs.R

inat <- read.csv("data/inat_phenology_data_subsp.csv", h=TRUE) %>% mutate(observed_on = ymd(observed_on))

glimpse(inat) # 11,082 raw observations
table(inat$phenology)
table(inat$phenology, inat$year)

# is fruit ever observed in Jan, Feb, or March? That's really the PREVIOUS flowering year
filter(inat, phenology=="Fruiting", month(observed_on)<4) # hmmm okay
inat <- filter(inat, !(phenology=="Fruiting" & month(observed_on)<4)) # clean those out


#-------------------------------------------------------------------------
# organize iNat observations for extraction of summarized PRISM data

# dealing with the 2019 "anomaly" by creating a second pseudo-year
inat$y2 <- inat$year
inat$y2[inat$year==2019 & month(inat$observed_on)>6] <- 2019.5

flowering <- data.frame(matrix(0,0,5))
names(flowering) <- c("lon","lat","year","type","flr")

prism_temp_rast <- raster("data/PRISM/annual/ppt_Mojave_2010Q1.bil")

# then ...
for(yr in 1:length(unique(inat$y2))){

yes <- 2*(rasterize(dplyr::filter(inat, y2==unique(inat$y2)[yr], phenology!="No Evidence of Flowering")[,c("longitude","latitude")], prism_temp_rast, fun=sum, background=0) > 0)

no <- rasterize(dplyr::filter(inat, y2==unique(inat$y2)[yr], phenology=="No Evidence of Flowering")[,c("longitude","latitude")], prism_temp_rast, fun=sum, background=0) > 0

yearyes <- rasterToPoints(yes+no, fun=function(x){x>1})
yearno <- rasterToPoints(yes+no, fun=function(x){x==1})

flowering <- rbind(flowering, data.frame(lon=c(yearyes[,"x"],yearno[,"x"]), lat=c(yearyes[,"y"],yearno[,"y"]), year=unique(inat$y2)[yr], flr=rep(c(TRUE,FALSE),c(nrow(yearyes),nrow(yearno)))))

}

head(flowering) # oh right we've lost type

# separate (sub)species again

# identify (sub)species
# inat_pheno_data <- read.csv("data/inat_phenology_data.csv", h=TRUE)
jtssps <- read_sf("data/Jotr_ssp_range.kml")

flo_in_ssp <- st_join(st_as_sf(flowering, coords=c("lon", "lat"), crs=crs(jtssps)), jtssps, join = st_within) %>% cbind(flowering[,c("lon", "lat")]) %>% as.data.frame(.) %>% dplyr::select(-geometry, -Description) %>% rename(type=Name) %>% dplyr::select(lon, lat, type, year, flr)

glimpse(flo_in_ssp) # okay okay okay!
table(flo_in_ssp$type, useNA="ifany")
table(flo_in_ssp$year, flo_in_ssp$flr) # think that looks good ...

write.table(flo_in_ssp, "output/flowering_obs_rasterized_subsp.csv", sep=",", col.names=TRUE, row.names=FALSE)

#-------------------------------------------------------------------------
# attach PRISM data to flowering/not flowering observations

# NEW predictor variables, informed by more specific hypotheses
# DEFINED: y0 is the year flowers are observed; y1 the year before, y2 two years before ...
# pptW0 - precip total, Oct y1 to Mar y0 --- "winter-of" precip
# pptY0 - precip total, Apr y1 to Mar y0 --- "year-of" precip
# pptW0W1 - precip total, Oct y2 to Mar y1 + Oct y1 to Mar y0 --- "winter accumulated" precip 
# pptY0W1 - precip total, Oct y2 to Mar y0 --- "year accumulated" precip 
# pptY0Y1 - precip total, Apr y2 to Mar y0 --- "all accumulated" precip 
# tmaxW0 - max temp, Oct y1 to Mar y0 --- "winter-of" max temp
# tminW0 - min temp, Oct y1 to Mar y0 --- "winter-of" min temp
# tmaxW0vW1 - difference in max temp, Oct y2 to Mar y1 vs Oct y1 to Mar y0 --- max-temp difference
# tminW0vW1 - difference in min temp, Oct y2 to Mar y1 vs Oct y1 to Mar y0 --- min-temp difference
# vpdmaxW0 - max VPD, Oct y1 to Mar y0 --- "winter-of" max VPD
# vpdminW0 - max VPD, Oct y1 to Mar y0 --- "winter-of" min VPD
# vpdmaxW0vW1 - difference in max VPD, Oct y2 to Mar y1 vs Oct y1 to Mar y0 --- max-VPD difference
# vpdminW0vW1 - difference in min VPD, Oct y2 to Mar y1 vs Oct y1 to Mar y0 --- min-VPD difference

flr.clim <- data.frame(matrix(0,0,ncol(flo_in_ssp)+13))
names(flr.clim) <- c(colnames(flo_in_ssp), "pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1", "vpdmaxW0", "vpdminW0", "vpdmaxW0vW1", "vpdminW0vW1")

# LOOP over years, because of the current year previous year thing ...
for(y in unique(flo_in_ssp$year)){

# y <- 2020 # test condition

# assemble weather data predictors for a given year
yr <- floor(y) # not going to try to lag 2019.5 for this predictor set 

preds <- brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_", yr,".grd", sep=""))

# pull subset of flowering observations for year
flsub <- subset(flo_in_ssp, year==y)
flsub <- cbind(flsub, raster::extract(preds, flsub[,c("lon","lat")], df=FALSE))

flr.clim <- rbind(flr.clim,flsub) 

write.table(flr.clim, "output/flowering_obs_climate_v2_subsp.csv", sep=",", col.names=TRUE, row.names=FALSE)

} # END LOOP over years


# and that's generated a data file we can feed into Embarcadero ... in the next script!








