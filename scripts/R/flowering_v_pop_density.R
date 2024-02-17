# Comparing predicted historical flowering to population density in Joshua tree
# Assumes local environment
# jby 2024.02.07

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
sdm.pres <- read_sf("../data/Yucca/Jotr_SDM2023_range.shp")
spp.ranges <- read_sf("data/Jotr_ssp_range.kml")

# Flowering years predicted by base and RI models, 
flyrs.jotr <- read.csv("output/jotr_reconstructed_flowering_years.csv")

flyrs_pred_change <- read.csv("output/jotr_flowering_predictors_change.csv") |> mutate(predictor = factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))


# Empirical maps
emp.yubr <- raster("../data/Yucca/joshua_tree_distribution/yubr_occupied_habitat.tif")
emp.yuja <- raster("../data/Yucca/joshua_tree_distribution/yuja_occupied_habitat.tif")

emp.jotr <- emp.yubr+emp.yuja # and now we have a two-species empirical layer

#-------------------------------------------------------------------------
# convert empirical map to density estimates

pred_rast <- stack("output/BART/jotr_BART_RI_predicted_flowering_1900-2023.gri") # raster to get the PRISM grid
projection(pred_rast) <- CRS("+init=epsg:4269")

emp.jotr.pts <- st_as_sf(data.frame(rasterToPoints(emp.jotr, fun=function(x){x==1})), coords=c("x", "y"), crs=crs(emp.jotr)) %>% st_transform(crs=crs(pred_rast))

jotr.density <- rasterize(emp.jotr.pts, field="layer", pred_rast, fun=sum)

pred.rast.density <- addLayer(pred_rast, jotr.density) # okay this confirms the whole shebang *should* be on a common grid and extent

jotr.density.df <- cbind(coordinates(jotr.density), as.data.frame(jotr.density)) %>% rename(lon=x, lat=y, pop_dens=layer) %>% filter(!is.na(pop_dens))

glimpse(jotr.density.df)
glimpse(flyrs.jotr)

#-------------------------------------------------------------------------
# merge and compare

flyrs.density <- flyrs.jotr %>% inner_join(jotr.density.df)
glimpse(flyrs.density) # 1,843 cells, but we expect there to be a bunch missing
filter(flyrs.density, is.na(pop_dens)) # check

cor.test(~pop_dens+flyrs_1900_1929, data=flyrs.density, method="sp", alt="greater") # n.s.
cor.test(~pop_dens+flyrs_1990_2019, data=flyrs.density, method="sp", alt="greater") # rho = 0.04, p = 0.06
cor.test(~pop_dens+flyrs_all, data=flyrs.density, method="sp", alt="greater") # n.s.

cor.test(~pop_dens+flyrs_change, data=flyrs.density, method="sp") # rho = 0.12, p = 1.017e-07


cor.test(~pop_dens+flyrs_RI_1900_1929, data=flyrs.density, method="sp", alt="greater") # n.s.
cor.test(~pop_dens+flyrs_RI_1990_2019, data=flyrs.density, method="sp", alt="greater") # n.s.
cor.test(~pop_dens+flyrs_RI_all, data=flyrs.density, method="sp", alt="greater") # n.s.

cor.test(~pop_dens+flyrs_RI_change, data=flyrs.density, method="sp") # n.s.


{cairo_pdf("output/figures/pop_density_vs_flyrs_change.pdf", width=4, height=4.5)

ggplot(flyrs.density, aes(y=flyrs_change, x=pop_dens/80)) + geom_point(alpha=0.5) + geom_smooth(method="lm") + 

labs(y = "Change in flowering years, 1990-2019 vs. 1900-1929", x = "Proportion area with Joshua trees present") + 
theme_bw()

}
dev.off()

hist(flyrs.density$flyrs_change)
hist(flyrs.density$pop_dens)
