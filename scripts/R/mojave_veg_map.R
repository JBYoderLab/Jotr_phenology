# Linking Central Mojave Vegetation Mapping Project plots to Joshua tree range
# Assumes local environment
# jby 2024.01.17

# starting up ------------------------------------------------------------

# setwd("~/Documents/Active_projects/Jotr_phenology")
# setwd("~/Documents/Academic/Active_projects/Jotr_phenology")

library("tidyverse")

library("raster")
library("sp")
library("sf")
library("hexbin")
library("embarcadero")

source("../shared/Rscripts/base.R") # my special mix of personal functions
source("../shared/Rscripts/base_graphics.R") # my special mix of personal functions


#-------------------------------------------------------------------------
# load up and organize data

# Jotr SDM and species boundaries
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range/Jotr_SDM2023_range.shp")
spp.ranges <- read_sf("data/Jotr_ssp_range.kml")

# MVMP plots as a point layer
mvmp_points <- read_sf("data/MVMP_plots_points/plots_points.shp")
glimpse(mvmp_points)

# MVMP "crosswalk" table
mvmp_cw <- read.csv("data/Mojave_Images_Crosswalk.csv")
glimpse(mvmp_cw)

mvmp_images <- filter(mvmp_cw, !grepl("^[Nn]o photos", PhotoID))
glimpse(mvmp_images) # that worked, good!

#-------------------------------------------------------------------------
# identify plot points within the JT range

ggplot() +
	geom_sf(data=sdm.pres) +
	geom_sf(data=mvmp_points)


mvmp_in_JTrange <- st_join(st_transform(mvmp_points, crs=crs(sdm.pres)), sdm.pres, join=st_within) %>% filter(!is.na(FID))

ggplot() +
	geom_sf(data=sdm.pres) +
	geom_sf(data=mvmp_points) +
	geom_sf(data=mvmp_in_JTrange, color="green")

mvmp_in_JTrange

#-------------------------------------------------------------------------
# identify sites with images in the JT range

mvmp_in_JTrange_photos <- mvmp_images %>% filter(FinalPlotCode %in% mvmp_in_JTrange$FinalPlotC)

glimpse(mvmp_in_JTrange_photos) # okay! 

# write out csv as a basis for our data collection sheet
write.table(mvmp_in_JTrange_photos[,1:4], "data/MEDP_survey_photos_JTrange.csv", sep=",", row.names=FALSE)


#-------------------------------------------------------------------------
# link coordinates with annotations and process for use as validation

annot <- read.csv("data/MDEP_annotations.csv")
glimpse(annot)
table(annot$FinalPlotCode) # lotta duplicates, eh?
table(annot$JT_present, (annot$Flowers | annot$Flower_budding | annot$Fruits)) # and very few flowers
  

annot_unique <- annot %>% filter(JT_present) %>% mutate(flr=(Flower_budding | Flowers | Fruits)) %>% group_by(FinalPlotCode, Survey_date) %>% summarize(flr=sum(flr)>0) %>% mutate(year=year(mdy(Survey_date)), month=month(mdy(Survey_date)))
glimpse(annot_unique)

mvmp_annot <- data.frame(st_coordinates(mvmp_in_JTrange), FinalPlotCode=mvmp_in_JTrange$FinalPlotC) %>% rename(lon=X, lat=Y) %>% right_join(annot_unique)

glimpse(mvmp_annot) # okay nice, about 143 more records?
table(mvmp_annot$flr, mvmp_annot$year)
table(mvmp_annot$month)

# adjust year of observation to flowering year
mvmp_annot$year[mvmp_annot$month<3] <- mvmp_annot$year[mvmp_annot$month<3]-1
table(mvmp_annot$month, mvmp_annot$year)

write.table(mvmp_annot, "output/MDEP_survey_records_cleaned.csv", sep=",")
