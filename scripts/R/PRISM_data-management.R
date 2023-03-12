# working with PRISM historical monthlys
# Assumes MAJEL environment 
# jby 2022.10.13

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")
library("lubridate")

library("raster")

library("prism")

if(!dir.exists("data/PRISM")) dir.create("data/PRISM")

prism_set_dl_dir("data/PRISM")

# Mojave crop extent (deliberately generous)
MojExt <- extent(-119, -112, 33, 38)


#-------------------------------------------------------------------------
# process PRISM data layers into a timeline archive I can use later

# make places to stash actual values and anomaly values
if(!dir.exists("data/PRISM/annual")) dir.create("data/PRISM/annual")


# FOR LOOP each year in PRISM
for(yr in 2019:2022){

# yr=2021

get_prism_monthlys(type="tmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="tmin", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="ppt", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmax", mon=1:12, year=yr, keepZip=FALSE)
get_prism_monthlys(type="vpdmin", mon=1:12, year=yr, keepZip=FALSE)


# FOR LOOP over quarters to do the thing ...
	for(q in 1:4){
	
	# q <- 1
	
	mos <- list(1:3,4:6,7:9,10:12)[[q]]
	
	# max temp in each quarter
	tmaxQ <- crop(max(raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmax", "monthly", year=yr, mon=mos[3])))), MojExt)
	
	writeRaster(tmaxQ, paste("data/PRISM/annual/tmax_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# min temp in each quarter
	tminQ <- crop(min(raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("tmin", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(tminQ, paste("data/PRISM/annual/tmin_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)
	
	# sum of precip in each quarter
	pptQ <- crop(sum(raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("ppt", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(pptQ, paste("data/PRISM/annual/ppt_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdmaxQ <- crop(sum(raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(vpdmaxQ, paste("data/PRISM/annual/vpdmax_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	vpdminQ <- crop(sum(raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(vpdminQ, paste("data/PRISM/annual/vpdmin_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	
	} # END loop over quarters

# clean out raw data
sapply(prism_archive_subset("tmax", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("tmin", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("ppt", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("vpdmax", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))
sapply(prism_archive_subset("vpdmin", "monthly", year=yr, mon=1:12), function(x) system(paste("rm -R data/PRISM/", x, sep="")))


} # END loop over years


#-------------------------------------------------------------------------
# normalize PRISM data to a 1981-2010 baseline

# make places to stash actual values and anomaly values
if(!dir.exists("data/PRISM/annual_normed_1981-2010")) dir.create("data/PRISM/annual_normed_1981-2010")
if(!dir.exists("data/PRISM/norms_1981-2010")) dir.create("data/PRISM/norms_1981-2010")

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


vpdmax_mn <- calc(brick(lapply(paste("data/PRISM/annual/vpdmax_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
vpdmax_sd <- calc(brick(lapply(paste("data/PRISM/annual/vpdmax_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(vpdmax_mn, paste("data/PRISM/norms_1981-2010/vpdmax_mean_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(vpdmax_sd, paste("data/PRISM/norms_1981-2010/vpdmax_sd_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)


vpdmin_mn <- calc(brick(lapply(paste("data/PRISM/annual/vpdmin_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), mean)
vpdmin_sd <- calc(brick(lapply(paste("data/PRISM/annual/vpdmin_Mojave_",1981:2010,"Q",q,".bil", sep=""), raster)), sd)

writeRaster(vpdmin_mn, paste("data/PRISM/norms_1981-2010/vpdmin_mean_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)
writeRaster(vpdmin_sd, paste("data/PRISM/norms_1981-2010/vpdmin_sd_Mojave_1981-2010_Q", q, ".bil", sep=""), overwrite=TRUE)

}


# FOR LOOP to norm the historical data ...
for(yr in 2019:2022){

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
	
	# max VPD in each quarter
	vpdmaxQnorm <- (raster(paste("data/PRISM/annual/vpdmax_Mojave_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/vpdmax_mean_Mojave_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/vpdmax_sd_Mojave_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(vpdmaxQnorm, paste("data/PRISM/annual_normed_1981-2010/vpdmax_normed_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)	

	# min VPD in each quarter
	vpdminQnorm <- (raster(paste("data/PRISM/annual/vpdmin_Mojave_", yr, "Q", q, ".bil", sep="")) - raster(paste("data/PRISM/norms_1981-2010/vpdmin_mean_Mojave_1981-2010_Q", q, ".bil", sep="")))/ raster(paste("data/PRISM/norms_1981-2010/vpdmin_sd_Mojave_1981-2010_Q", q, ".bil", sep=""))
	
	writeRaster(vpdminQnorm, paste("data/PRISM/annual_normed_1981-2010/vpdmin_normed_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)	

	
	} # END loop over quarters

}

# and now I have quarterly historical data for the Mojave to work with for all downstream analysis


#-------------------------------------------------------------------------
# calculate hypothesized predictors from normalized PRISM data

# make places to stash actual values and anomaly values
if(!dir.exists("data/PRISM/derived_predictors")) dir.create("data/PRISM/derived_predictors")

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


# FOR LOOP over the historical data ...
# nb, because retrospective, can't do all the way back to 1895
for(yr in 1900:2022){

# yr <- 1900

pptW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep="")), raster)),sum)
pptY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:2,".bil", sep="")), raster)),sum)
pptW0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q4.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0W1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q4.bil", sep="")), raster)),sum)
pptY0Y1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:1,".bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q",4:2,".bil", sep="")), raster)),sum)
tmaxW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q4.bil", sep="")), raster)), max)
tminW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q4.bil", sep="")), raster)), min)
tmaxW0vW1 <- tmaxW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-2,"Q4.bil", sep="")), raster)), max)
tminW0vW1 <- tminW0 - calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-2,"Q4.bil", sep="")), raster)), min)
vpdmaxW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmax_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmax_Mojave_",yr-1,"Q4.bil", sep="")), raster)), max)
vpdminW0 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmin_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmin_Mojave_",yr-1,"Q4.bil", sep="")), raster)), min)
vpdmaxW0vW1 <- vpdmaxW0 - calc(brick(lapply(c(paste("data/PRISM/annual/vpdmax_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmax_Mojave_",yr-2,"Q4.bil", sep="")), raster)), max)
vpdminW0vW1 <- vpdminW0 - calc(brick(lapply(c(paste("data/PRISM/annual/vpdmin_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmin_Mojave_",yr-2,"Q4.bil", sep="")), raster)), min)

preds <- brick(c(pptW0, pptY0, pptW0W1, pptY0W1, pptY0Y1, tmaxW0, tminW0, tmaxW0vW1, tminW0vW1, vpdmaxW0, vpdminW0, vpdmaxW0vW1, vpdminW0vW1))
names(preds) <- c("pptW0", "pptY0", "pptW0W1", "pptY0W1", "pptY0Y1", "tmaxW0", "tminW0", "tmaxW0vW1", "tminW0vW1", "vpdmaxW0", "vpdminW0", "vpdmaxW0vW1", "vpdminW0vW1")

writeRaster(preds, paste("data/PRISM/derived_predictors/PRISM_derived_predictors_", yr,".grd", sep=""), overwrite=TRUE) # confirmed write-out and read-in

} # and now each year has the predictor set stashed as a multilayer raster, cool

