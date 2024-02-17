# working with PRISM historical monthlys
# Assumes MAJEL environment 
# jby 2024.02.16

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
MojExt <- extent(-119, -112, 33, 40)


#-------------------------------------------------------------------------
# process PRISM data layers into a timeline archive I can use later

# make places to stash actual values and anomaly values
if(!dir.exists("data/PRISM/annual")) dir.create("data/PRISM/annual")


# FOR LOOP each year in PRISM
for(yr in 1895:2021){

# yr=2023

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

	# max vpd in each quarter
	vpdmaxQ <- crop(max(raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmax", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(vpdmaxQ, paste("data/PRISM/annual/vpdmax_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	# min vpd in each quarter
	vpdminQ <- crop(min(raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[1]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[2]))), raster(pd_to_file(prism_archive_subset("vpdmin", "monthly", year=yr, mon=mos[3])))), MojExt)

	writeRaster(vpdminQ, paste("data/PRISM/annual/vpdmin_Mojave_", yr, "Q", q, ".bil", sep=""), overwrite=TRUE)

	# re-insert file clearance here!!
	
	} # END loop over quarters


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
# DEFINED: y0 is Apr-Mar of the year flowers are observed; y1 the year before, y2 two years before ...
# for each of YEAR 0, 1, and 2 ...
# total precip (ppt)
# max and min temperature (tmax and tmin)
# max and min vapor pressure deficit (vpdmax and vpdmin)
# and then also CONTRASTS Y0-Y1, Y1-Y2 for each of these

# FOR LOOP over the historical data ...
# nb, because retrospective, can't do all the way back to 1895
rast_ref <- raster("data/PRISM/annual/ppt_Mojave_2015Q1.bil") # standard to resample against, because WTF

for(yr in 1900:2022){

# yr <- 2023

# ppt, all summarized
pptY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),sum) # borked in 2021
pptY1 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),sum)
pptY2 <- calc(brick(lapply(c(paste("data/PRISM/annual/ppt_Mojave_",yr-2,"Q1.bil", sep=""), paste("data/PRISM/annual/ppt_Mojave_",yr-3,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),sum)

pptY0Y1 <- pptY0 - pptY1
pptY1Y2 <- pptY1 - pptY2

# tmax, all summarized
tmaxY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),max) # borked in 2021
tmaxY1 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-2,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),max)
tmaxY2 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmax_Mojave_",yr-2,"Q1.bil", sep=""), paste("data/PRISM/annual/tmax_Mojave_",yr-3,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),max)

tmaxY0Y1 <- tmaxY0 - tmaxY1
tmaxY1Y2 <- tmaxY1 - tmaxY2

# tmin, all summarized
tminY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),min) # borked in 2021
tminY1 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-2,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),min)
tminY2 <- calc(brick(lapply(c(paste("data/PRISM/annual/tmin_Mojave_",yr-2,"Q1.bil", sep=""), paste("data/PRISM/annual/tmin_Mojave_",yr-3,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),min)

tminY0Y1 <- tminY0 - tminY1
tminY1Y2 <- tminY1 - tminY2

# vpdmax, all summarized
vpdmaxY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmax_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmax_Mojave_",yr-1,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),max) # borked in 2021
vpdmaxY1 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmax_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmax_Mojave_",yr-2,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),max)
vpdmaxY2 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmax_Mojave_",yr-2,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmax_Mojave_",yr-3,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),max)

vpdmaxY0Y1 <- vpdmaxY0 - vpdmaxY1
vpdmaxY1Y2 <- vpdmaxY1 - vpdmaxY2

# vpdmin, all summarized
vpdminY0 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmin_Mojave_",yr,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmin_Mojave_",yr-1,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),min) # borked in 2021
vpdminY1 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmin_Mojave_",yr-1,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmin_Mojave_",yr-2,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),min)
vpdminY2 <- calc(brick(lapply(c(paste("data/PRISM/annual/vpdmin_Mojave_",yr-2,"Q1.bil", sep=""), paste("data/PRISM/annual/vpdmin_Mojave_",yr-3,"Q",4:2,".bil", sep="")), function(x) resample(raster(x), rast_ref))),min)

vpdminY0Y1 <- vpdminY0 - vpdminY1
vpdminY1Y2 <- vpdminY1 - vpdminY2

preds <- brick(c(pptY0, pptY1, pptY2, pptY0Y1, pptY1Y2, tmaxY0, tmaxY1, tmaxY2, tmaxY0Y1, tmaxY1Y2, tminY0, tminY1, tminY2, tminY0Y1, tminY1Y2, vpdmaxY0, vpdmaxY1, vpdmaxY2, vpdmaxY0Y1, vpdmaxY1Y2, vpdminY0, vpdminY1, vpdminY2, vpdminY0Y1, vpdminY1Y2))
names(preds) <- c("pptY0", "pptY1", "pptY2", "pptY0Y1", "pptY1Y2", "tmaxY0", "tmaxY1", "tmaxY2", "tmaxY0Y1", "tmaxY1Y2", "tminY0", "tminY1", "tminY2", "tminY0Y1", "tminY1Y2", "vpdmaxY0", "vpdmaxY1", "vpdmaxY2", "vpdmaxY0Y1", "vpdmaxY1Y2", "vpdminY0", "vpdminY1", "vpdminY2", "vpdminY0Y1", "vpdminY1Y2")

writeRaster(preds, paste("data/PRISM/derived_predictors/PRISM_derived_predictors_", yr,".grd", sep=""), overwrite=TRUE) # confirmed write-out and read-in

} # and now each year has the predictor set stashed as a multilayer raster, cool


#-------------------------------------------------------------------------
# summarize predictors, for SI

library("sf")

sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range.shp")
sdm.buff <- st_buffer(st_transform(sdm.pres[,2], crs=3857), 1000) # put a 1km buffer on the range polygons


jotr.maskHist <- mask(jotr.histStack, st_transform(sdm.buff, crs=4269), touches=TRUE)

# this may take a moment ...
predsComp <- rbind(
		cbind(year=2008, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2008.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2009, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2009.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2010, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2010.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2011, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2011.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2012, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2012.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2013, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2013.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2014, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2014.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2015, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2015.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2016, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2016.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2017, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2017.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2018, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2018.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2019, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2019.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2020, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2020.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2021, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2021.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0))),
		cbind(year=2022, as.data.frame(mask(stack("data/PRISM/derived_predictors/PRISM_derived_predictors_2022.grd"), st_transform(sdm.buff, crs=4269), touches=TRUE)) %>% filter(!is.na(pptY0)))
		)

write.table(predsComp, "output/derived_predictors_jtrange.csv", sep=",", col.names=TRUE, row.names=FALSE)

preds.ln <- predsComp %>% pivot_longer(all_of(c("pptY0", "pptY1", "pptY2", "pptY0Y1", "pptY1Y2", "tmaxY0", "tmaxY1", "tmaxY2", "tmaxY0Y1", "tmaxY1Y2", "tminY0", "tminY1", "tminY2", "tminY0Y1", "tminY1Y2", "vpdmaxY0", "vpdmaxY1", "vpdmaxY2", "vpdmaxY0Y1", "vpdmaxY1Y2", "vpdminY0", "vpdminY1", "vpdminY2", "vpdminY0Y1", "vpdminY1Y2")), names_to="predictor", values_to="value")

glimpse(preds.ln)

# now a summary figure?
xnames <- c("pptY0", "pptY1", "pptY2", "pptY0Y1", "pptY1Y2", "tmaxY0", "tmaxY0Y1", "tminY0", "tminY0Y1", "vpdmaxY0", "vpdmaxY0Y1", "vpdminY0", "vpdminY0Y1")

# c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

# or a figure?
jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

preds.fig <- filter(preds.ln, predictor%in%jotr.preds)
preds.fig$predictor[preds.fig$predictor=="pptY1Y2"] <- "Delta[Y1-2]*PPT"
preds.fig$predictor[preds.fig$predictor=="pptY0Y1"] <- "Delta[Y0-1]*PPT"
preds.fig$predictor[preds.fig$predictor=="vpdmaxY0"] <- "Max*VPD[Y0]"
preds.fig$predictor[preds.fig$predictor=="vpdminY0Y1"] <- "Delta[Y0-1]*Min*VPD"
preds.fig$predictor[preds.fig$predictor=="tminY0"] <- "Min*Temp[Y0]"
preds.fig$predictor[preds.fig$predictor=="tmaxY0Y1"] <- "Delta[Y0-1]*Max*Temp"

preds.fig$predictor <- as.factor(preds.fig$predictor, c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp"))

{cairo_pdf("output/figures/SI_predictors.pdf", width=6.5, height=6)

ggplot(preds.fig, aes(x=year, y=value, group=year)) + geom_boxplot(outlier.alpha=0.5, outlier.size=0.5) + 

facet_wrap("predictor", scale="free_y", ncol=2, labeller="label_parsed") +

labs(y="Predictor value", x="Year") +

theme_bw()

}
dev.off()

# maybe a summary table?
preds.ln %>% filter(predictor %in% xnames) %>% group_by(predictor, year) %>% summarize(mn.annual = mean(value)) %>% group_by(predictor) %>% summarize(mnPred = mean(mn.annual), st.error = sd(mn.annual)/sqrt(15))


# and in training data specifically

flow <- read.csv("output/flowering_obs_climate_subsp.csv") # flowering/not flowering, biologically-informed candidate predictors, subspecies id'd

dim(flow)
glimpse(flow)


# variant datasets -- dealing with the second flowering in 2018
flow2 <- flow %>% filter(!(year==2018.5 & flr==TRUE), year>=2008, year<2023) %>% mutate(year=floor(year)) # drop the late-flowering anomaly
glimpse(flow2) # 2,632 in our final working set

flow2 %>% dplyr::select(lat, lon, year, all_of(xnames)) %>% pivot_longer(all_of(xnames), names_to="predictor", values_to="value") %>% group_by(predictor, year) %>% summarize(mn.annual = mean(value)) %>% group_by(predictor) %>% summarize(mnPred = mean(mn.annual), st.error = sd(mn.annual)/sqrt(15))
