# working with phenology-annotated iNat observations
# Assumes MAJEL environment 
# jby 2022.07.05

# starting up ------------------------------------------------------------

# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("lubridate")

library("raster")

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
# start with monthly tmin, tmax, and ppt for the year of flowering and the year previous

flr.clim <- data.frame(matrix(0,0,ncol(flowering)+4*3*3))
names(flr.clim) <- c(colnames(flowering), paste("ppt",paste("y",rep(0,4),sep=""),"q",1:4,sep=""), paste("tmax",paste("y",rep(0,4),sep=""),"q",1:4,sep=""), paste("tmin",paste("y",rep(0,4),sep=""),"q",1:4,sep=""), paste("ppt",paste("y",rep(1,4),sep=""),"q",1:4,sep=""), paste("tmax",paste("y",rep(1,4),sep=""),"q",1:4,sep=""), paste("tmin",paste("y",rep(1,4),sep=""),"q",1:4,sep=""), paste("ppt",paste("y",rep(2,4),sep=""),"q",1:4,sep=""), paste("tmax",paste("y",rep(2,4),sep=""),"q",1:4,sep=""), paste("tmin",paste("y",rep(2,4),sep=""),"q",1:4,sep=""))

# LOOP over years, because of the current year previous year thing ...
for(yr in unique(flowering$year)){

# for each year's observations, pull monthly climate values from y0 and y-1
# yr = 2020

flsub <- subset(flowering, year==yr)

if(yr!=2019.5) files <- c(list.files("data/PRISM/annual", pattern=paste(yr,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual", pattern=paste(yr-1,"Q\\d","\\.bil",sep=""), full=TRUE),  list.files("data/PRISM/annual", pattern=paste(yr-2,"Q\\d","\\.bil",sep=""), full=TRUE))

if(yr==2019.5) files <- c(list.files("data/PRISM/annual", pattern=paste(yr-0.5,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual", pattern=paste(yr-1.5,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual", pattern=paste(yr-2.5,"Q\\d","\\.bil",sep=""), full=TRUE), list.files("data/PRISM/annual", pattern=paste(yr-3.5,"Q\\d","\\.bil",sep=""), full=TRUE))[c(3,4,13,14, 7,8,17,18, 11,12,21,22, 15,16,25,26, 19,20,29,30, 23,24,33,34, 27,28,37,38, 31,32,41,42, 35,36,45,46)]


# LOOP over the files for data extraction
for(f in files){

# f = files0[1]

flsub <- cbind(flsub, raster::extract(raster(f), flsub[,c("lon","lat")], df=FALSE))

} # END LOOP over files

colnames(flsub) <- colnames(flr.clim)

flr.clim <- rbind(flr.clim,flsub)

write.table(flr.clim, "output/flowering_obs_climate.csv", sep=",", col.names=TRUE, row.names=FALSE)

} # END LOOP over years


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





