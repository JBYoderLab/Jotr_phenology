# working with PRISM historical monthlys
# Assumes MAJEL environment 
# jby 2022.07.12

# starting up ------------------------------------------------------------

# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("lubridate")

library("raster")

library("prism")

prism_set_dl_dir("data/PRISM")

if(!dir.exists("data/PRISM")) dir.create("data/PRISM")

# Mojave crop extent (deliberately generous)
MojExt <- extent(-119, -112, 33, 38)

#-------------------------------------------------------------------------
# process PRISM normals

# first, downloads
get_prism_normals(type="tmax", resolution="4km", mon=1:12, keepZip=FALSE)
get_prism_normals(type="tmin", resolution="4km", mon=1:12, keepZip=FALSE)
get_prism_normals(type="ppt", resolution="4km", mon=1:12, keepZip=FALSE)

# next, make a place to stash them
if(!dir.exists("data/PRISM/normals")) dir.create("data/PRISM/normals")

# identify the file sets I'm working with
tmaxNorms <- grep("tmax_30yr_normal", prism_archive_ls(), value=TRUE)
tminNorms <- grep("tmin_30yr_normal", prism_archive_ls(), value=TRUE)
pptNorms <- grep("ppt_30yr_normal", prism_archive_ls(), value=TRUE)

# loop over quarters to do the thing ...
for(q in 1:4){

# q <- 1

mos <- list(1:3,4:6,7:9,10:12)[[q]]

# max temp in each quarter
writeRaster(crop(max(raster(paste("data/PRISM/", tmaxNorms[mos[1]], "/", tmaxNorms[mos[1]], ".bil", sep="")), raster(paste("data/PRISM/", tmaxNorms[mos[2]], "/", tmaxNorms[mos[2]], ".bil", sep="")), raster(paste("data/PRISM/", tmaxNorms[mos[3]], "/", tmaxNorms[mos[3]], ".bil", sep=""))), MojExt), paste("data/PRISM/normals/tmax_Mojave_30yr_normal_4km_Q", q, ".bil", sep=""), overwrite=TRUE)

# min temp in each quarter
writeRaster(crop(min(raster(paste("data/PRISM/", tminNorms[mos[1]], "/", tminNorms[mos[1]], ".bil", sep="")), raster(paste("data/PRISM/", tminNorms[mos[2]], "/", tminNorms[mos[2]], ".bil", sep="")), raster(paste("data/PRISM/", tminNorms[mos[3]], "/", tminNorms[mos[3]], ".bil", sep=""))), MojExt), paste("data/PRISM/normals/tmin_Mojave_30yr_normal_4km_Q", q, ".bil", sep=""), overwrite=TRUE)

# sum of precip in each quarter
writeRaster(crop(sum(raster(paste("data/PRISM/", pptNorms[mos[1]], "/", pptNorms[mos[1]], ".bil", sep="")), raster(paste("data/PRISM/", pptNorms[mos[2]], "/", pptNorms[mos[2]], ".bil", sep="")), raster(paste("data/PRISM/", pptNorms[mos[3]], "/", pptNorms[mos[3]], ".bil", sep=""))), MojExt), paste("data/PRISM/normals/ppt_Mojave_30yr_normal_4km_Q", q, ".bil", sep=""), overwrite=TRUE)

}


#-------------------------------------------------------------------------
# process PRISM normals into a timeline archive I can use later

# make places to stash actual values and anomaly values
if(!dir.exists("data/PRISM/annual")) dir.create("data/PRISM/annual")


# FOR LOOP each year in PRISM
for(yr in 1895:2022){

# yr=2021

get_prism_monthlys(type="tmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmin", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="ppt", mon=1:12, year=yr, keepZip=FALSE)


# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){
	
	# q <- 1
	
	mos <- list(1:3,4:6,7:9,10:12)[[q]]
	
	# 1981-2010 normals rasters
	tmaxNorm <- raster(paste("data/PRISM/normals/tmax_Mojave_30yr_normal_4km_Q", q, ".bil", sep=""))
	tminNorm <- raster(paste("data/PRISM/normals/tmin_Mojave_30yr_normal_4km_Q", q, ".bil", sep=""))
	pptNorm <- raster(paste("data/PRISM/normals/ppt_Mojave_30yr_normal_4km_Q", q, ".bil", sep=""))
	
	# max temp in each quarter
	tmaxQ <- crop(max(raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[3])))), MojExt)
	
	writeRaster(tmaxQ, paste("data/PRISM/annual/tmax_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# min temp in each quarter
	tminQ <- crop(min(raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(tminQ, paste("data/PRISM/annual/tmin_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# sum of precip in each quarter
	pptQ <- crop(sum(raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(pptQ, paste("data/PRISM/annual/ppt_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	} # END loop over quarters

# clean out raw data
sapply(prism_archive_subset("tmax", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("tmin", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("ppt", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))

} # END loop over years


#-------------------------------------------------------------------------
# process PRISM data into a timeline archive I can use later
# but this time normalized to a 1981-2010 baseline

# make places to stash actual values and anomaly values
dir.create("data/PRISM/annual_normed_1981-2010")
dir.create("data/PRISM/norms_1981-2010")

quarterlies <- list.files("data/PRISM/annual", pattern=".bil")

# first, our normals, 1981-2010
for(q in 1:4){

# q = 1

tmax_mn <- calc(brick(lapply(paste("data/PRISM/annual/tmax_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
tmax_sd <- calc(brick(lapply(paste("data/PRISM/annual/tmax_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(tmax_mn, paste("data/PRISM/norms_1981-2010/tmax_mean_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(tmax_sd, paste("data/PRISM/norms_1981-2010/tmax_sd_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)


tmin_mn <- calc(brick(lapply(paste("data/PRISM/annual/tmin_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
tmin_sd <- calc(brick(lapply(paste("data/PRISM/annual/tmin_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(tmin_mn, paste("data/PRISM/norms_1981-2010/tmin_mean_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(tmin_sd, paste("data/PRISM/norms_1981-2010/tmin_sd_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)


ppt_mn <- calc(brick(lapply(paste("data/PRISM/annual/ppt_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
ppt_sd <- calc(brick(lapply(paste("data/PRISM/annual/ppt_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(ppt_mn, paste("data/PRISM/norms_1981-2010/ppt_mean_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(ppt_sd, paste("data/PRISM/norms_1981-2010/ppt_sd_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)

}


# FOR LOOP to norm the historical data ...
for(yr in 1895:2022){

# yr=1895

# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){
	
	# q <- 1
	
	# max temp in each quarter
	tmaxQnorm <- (raster(paste("data/PRISM/annual/tmax_Mojave_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/tmax_mean_Mojave_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/tmax_sd_Mojave_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(tmaxQnorm, paste("data/PRISM/annual_normed_1981-2010/tmax_normed_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# min temp in each quarter
	tminQnorm <- (raster(paste("data/PRISM/annual/tmin_Mojave_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/tmin_mean_Mojave_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/tmin_sd_Mojave_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(tminQnorm, paste("data/PRISM/annual_normed_1981-2010/tmin_normed_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# sum of precip in each quarter
	pptQnorm <- crop((raster(paste("data/PRISM/annual/ppt_Mojave_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/ppt_mean_Mojave_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/ppt_sd_Mojave_1981-2010_Q", q, ".bil", sep="")), MojExt)
	
	writeRaster(pptQnorm, paste("data/PRISM/annual_normed_1981-2010/ppt_normed_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	} # END loop over quarters

}


# and now I have quarterly historical data for the Mojave to work with for all downstream analysis

# need to re-run this for 1895 tomorrow



