# working with phenology-annotated iNat observations
# Assumes MAJEL environment 
# jby 2022.07.12

# starting up ------------------------------------------------------------

# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("lubridate")

library("raster")
library("sf")

# Mojave crop extent (deliberately generous)
MojExt <- extent(-119, -112, 33, 38)

#-------------------------------------------------------------------------
# read in iNat observations compiled using inat_phenology_obs.R

inat <- read.csv("data/inat_phenology_data.csv", h=TRUE) %>% mutate(observed_on = ymd(observed_on))

glimpse(inat)
table(inat$phenology)
table(inat$phenology, inat$year)

# is fruit ever observed in Jan, Feb, or March? That's really the PREVIOUS flowering year
inat$year[inat$phenology=="Fruiting" && month(inat$observed_on)<4] # okay then


#-------------------------------------------------------------------------
# organize iNat observations for extraction of summarized PRISM data

# dealing with the 2019 "anomaly" by creating a second pseudo-year
inat$y2 <- inat$year
inat$y2[inat$year==2019 & month(inat$observed_on)>6] <- 2019.5

flowering <- data.frame(matrix(0,0,4))
names(flowering) <- c("lon","lat","year","flr")

prism_temp_rast <- raster("data/PRISM/annual/ppt_Mojave_2010Q1.bil")

# then ...
for(yr in 1:length(unique(inat$y2))){

yes <- 2*(rasterize(dplyr::filter(inat, y2==unique(inat$y2)[yr], phenology!="No Evidence of Flowering")[,c("longitude","latitude")], prism_temp_rast, fun=sum, background=0) > 0)

no <- rasterize(dplyr::filter(inat, y2==unique(inat$y2)[yr], phenology=="No Evidence of Flowering")[,c("longitude","latitude")], prism_temp_rast, fun=sum, background=0) > 0

yearyes <- rasterToPoints(yes+no, fun=function(x){x>1})
yearno <- rasterToPoints(yes+no, fun=function(x){x==1})

flowering <- rbind(flowering, data.frame(lon=c(yearyes[,"x"],yearno[,"x"]), lat=c(yearyes[,"y"],yearno[,"y"]), year=unique(inat$y2)[yr], flr=rep(c(TRUE,FALSE),c(nrow(yearyes),nrow(yearno)))))

}

head(flowering)
glimpse(flowering) # okay okay okay!
table(flowering$year, flowering$flr) # think that looks good ...

write.table(flowering, "output/flowering_obs_rasterized.csv", sep=",", col.names=TRUE, row.names=FALSE)

#-------------------------------------------------------------------------
# attach PRISM data to flowering/not flowering observations

# we're going to generate some wiiiiiiide tables here
# NEW predictor variables, informed by more specific hypotheses
# DEFINED: y0 is the year flowers are observed; y1 the year before, y2 two years before ...
# - precip total, Oct y1 to Mar y0 --- "winter-of" precip (covers winter before flowering)
# - precip total, Apr y1 to Mar y0 --- "year-of" precip (covers summer + winter before)
# - precip total, Oct y2 to Mar y1 + Oct y1 to Mar y0 --- "winter accumulated" precip (covers two winters)
# - precip total, Oct y2 to Mar y0 --- "year accumulated" precip (covers two full winters and the summer between)
# - precip total, Apr y2 to Mar y0 --- "all accumulated" precip (covers two full summers and winters)
# - max temp, Oct y1 to Mar y0 --- "winter-of" max temp
# - min temp, Oct y1 to Mar y0 --- "winter-of" min temp
# - difference in max temp, Oct y2 to Mar y1 vs Oct y1 to Mar y0 --- max-temp difference
# - difference in min temp, Oct y2 to Mar y1 vs Oct y1 to Mar y0 --- min-temp difference

flr.clim <- data.frame(matrix(0,0,ncol(flowering)+9))
names(flr.clim) <- c(colnames(flowering), "pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1")

# LOOP over years, because of the current year previous year thing ...
for(y in unique(flowering$year)){

# y <- 2020 # test condition

# assemble weather data predictors for a given year
yr <- floor(y) # not going to try to lag 2019.5 for this predictor set 

pptW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep="")), raster)),sum)
pptY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:2,".bil", sep="")), raster)),sum)
pptW0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0Y1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q",4:2,".bil", sep="")), raster)),sum)
tmaxW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q4.bil", sep="")), raster)), max)
tminW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q4.bil", sep="")), raster)), min)
tmaxW0vW1 <- tmaxW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-2,"Q4.bil", sep="")), raster)), max)
tminW0vW1 <- tminW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-2,"Q4.bil", sep="")), raster)), min)

preds <- brick(c(pptW0, pptY0, pptW0W1, pptY0W1, pptY0Y1, tmaxW0, tminW0, tmaxW0vW1, tminW0vW1))
names(preds) <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1")

# pull subset of flowering observations for year
flsub <- subset(flowering, year==y)
flsub <- cbind(flsub, raster::extract(preds, flsub[,c("lon","lat")], df=FALSE))

flr.clim <- rbind(flr.clim,flsub) 

write.table(flr.clim, "output/flowering_obs_climate_v2.csv", sep=",", col.names=TRUE, row.names=FALSE)

} # END LOOP over years

# identify (sub)species
# flr.clim <- read.csv("output/flowering_obs_climate_v2.csv", h=TRUE)
jtssps <- read_sf("data/Jotr_range.kml")

obs_in_ssp <- st_join(st_as_sf(flr.clim, coords=c("lon", "lat"), crs=crs(jtssps)), jtssps, join = st_within) %>% cbind(flr.clim[,c("lon","lat")]) %>% as.data.frame(.) %>% dplyr::select(-geometry, -Description) %>% mutate(Name = gsub("(\\w+) Joshua tree range", "\\1", Name)) %>% rename(type=Name) %>% dplyr::select(lat, lon, type, year:tminW0vW1)

write.table(obs_in_ssp, "output/flowering_obs_climate_v2_subsp.csv", sep=",", col.names=TRUE, row.names=FALSE)

# and that's generated a data file we can feed into Embarcadero ... in a new script!



# same thing, but with values normalized to 1981-2010

flr.clim <- data.frame(matrix(0,0,ncol(flowering)+4*3*3))
names(flr.clim) <- c(colnames(flowering), paste("ppt",paste("y",rep(0,4),sep=""),"q",1:4,sep=""), paste("tmax",paste("y",rep(0,4),sep=""),"q",1:4,sep=""), paste("tmin",paste("y",rep(0,4),sep=""),"q",1:4,sep=""), paste("ppt",paste("y",rep(1,4),sep=""),"q",1:4,sep=""), paste("tmax",paste("y",rep(1,4),sep=""),"q",1:4,sep=""), paste("tmin",paste("y",rep(1,4),sep=""),"q",1:4,sep=""), paste("ppt",paste("y",rep(2,4),sep=""),"q",1:4,sep=""), paste("tmax",paste("y",rep(2,4),sep=""),"q",1:4,sep=""), paste("tmin",paste("y",rep(2,4),sep=""),"q",1:4,sep=""))

# LOOP over years, because of the current year previous year thing ...
for(yr in unique(flowering$year)){

# for each year's observations, pull monthly climate values from y0 and y-1
# yr = 2020

flsub <- subset(flowering, year==yr)

if(yr!=2019.5) files <- c(list.files("data/PRISM/annual_normed_1981-2010", pattern=paste(yr,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual_normed_1981-2010", pattern=paste(yr-1,"Q\\d","\\.bil",sep=""), full=TRUE),  list.files("data/PRISM/annual_normed_1981-2010", pattern=paste(yr-2,"Q\\d","\\.bil",sep=""), full=TRUE))

if(yr==2019.5) files <- c(list.files("data/PRISM/annual_normed_1981-2010", pattern=paste(yr-0.5,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual_normed_1981-2010", pattern=paste(yr-1.5,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual_normed_1981-2010", pattern=paste(yr-2.5,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual_normed_1981-2010", pattern=paste(yr-3.5,"Q\\d","\\.bil",sep=""), full=TRUE))[c(3,4,13,14, 7,8,17,18, 11,12,21,22, 15,16,25,26, 19,20,29,30, 23,24,33,34, 27,28,37,38, 31,32,41,42, 35,36,45,46)]


# LOOP over the files for data extraction
for(f in files){

# f = files0[1]

flsub <- cbind(flsub, raster::extract(raster(f), flsub[,c("lon","lat")], df=FALSE))

} # END LOOP over files

colnames(flsub) <- colnames(flr.clim)

flr.clim <- rbind(flr.clim,flsub)

write.table(flr.clim, "output/flowering_obs_climate_normed.csv", sep=",", col.names=TRUE, row.names=FALSE)

} # END LOOP over years





