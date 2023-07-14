# Analyzing predicted historical flowering in Joshua tree
# Assumes local environment
# jby 2023.07.05

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

library("ggdark")

#-------------------------------------------------------------------------
# load up and organize historical data

# Jotr SDM and species boundaries
sdm.pres <- read_sf("../data/Yucca/jotr_BART_sdm_pres", "jotr_BART_sdm_pres")
spp.ranges <- read_sf("data/Jotr_ssp_range.kml")

yubr.pres <- st_zm(st_intersection(sdm.pres, st_transform(filter(spp.ranges, Name=="YUBR"), crs=crs(sdm.pres))))
yuja.pres <- st_zm(st_intersection(sdm.pres, st_transform(filter(spp.ranges, Name=="YUJA"), crs=crs(sdm.pres))))

# flowering observation data
obs <- read.csv("output/flowering_obs_climate_subsp.csv") %>% mutate(y2 = year) %>% mutate(year=floor(year), type = factor(type, c("YUBR", "YUJA")))

# raster files of predicted prFL
jotr.files <- list.files("output/BART/predictions", pattern=".bil", full=TRUE)
jotr.ri.files <- list.files("output/BART/RI.predictions", pattern=".bil", full=TRUE)

yubr.files <- list.files("output/BART/predictions.YUBR", pattern=".bil", full=TRUE)
yuja.files <- list.files("output/BART/predictions.YUJA", pattern=".bil", full=TRUE)

MojExt <- extent(-119, -112, 33, 38) # Mojave extent, maybe useful

# expert validation observations
vobs <- read.csv("data/Validation_obs_by_spp.csv")

glimpse(vobs)

table(vobs$obs_by)

jotr.preds <- c("pptY1Y2", "pptY0Y1", "vpdmaxY0", "vpdminY0Y1", "tminY0", "tmaxY0Y1")

# Jotr, without RI --------------------------------
jotr.histStack <- raster::stack(sapply(jotr.files, function(x) crop(raster::raster(x), MojExt)))
names(jotr.histStack) <- paste("prFL",1900:2022,sep=".")
projection(jotr.histStack)<-CRS("+init=epsg:4269")

jotr.histStack

jotr.maskHist <- mask(jotr.histStack, st_transform(sdm.pres[,2], crs=4269))

jotr.hist.flowering <- cbind(coordinates(jotr.maskHist), as.data.frame(jotr.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(jotr.hist.flowering)

# Jotr, with RI ------------------------------------
jotr.ri.histStack <- raster::stack(sapply(jotr.ri.files, function(x) crop(raster::raster(x), MojExt)))
names(jotr.ri.histStack) <- paste("prFL",1900:2022,sep=".")
projection(jotr.ri.histStack)<-CRS("+init=epsg:4269")

jotr.ri.histStack

jotr.ri.maskHist <- mask(jotr.ri.histStack, st_transform(sdm.pres[,2], crs=4269))

jotr.ri.hist.flowering <- cbind(coordinates(jotr.ri.maskHist), as.data.frame(jotr.ri.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="RI.prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(jotr.ri.hist.flowering)

jotr.flowering.histories <- jotr.hist.flowering |> left_join(jotr.ri.hist.flowering)

glimpse(jotr.flowering.histories)

write.table(jotr.flowering.histories, "output/historic_flowering_reconst_jotr.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE) # write out and read in

jotr.hist.flowering <- read.csv("output/historic_flowering_reconst_jotr.csv")

# comparing predictions with and without RI
RItest <- jotr.hist.flowering |> filter(year>=2008)

baseVri <- table(RItest$prFL>0.26, RItest$RI.prFL>0.25) 
baseVri # rows, prFL; columns, RI.prFL

sum(baseVri[c(1,4)])/sum(baseVri) # overall 0.84 agreed
baseVri[2]/sum(baseVri) # base says T, RI says F
baseVri[3]/sum(baseVri) # RI says T, base says F

chisq.test(table(RItest$prFL>0.26, RItest$RI.prFL>0.25)) # p < 2.2e-16, natch

RIvNon <- RItest |> group_by(year) |> do(broom::tidy(cor.test(~prFL+RI.prFL, data=., method="spearman"))) %>% ungroup()

glimpse(RIvNon)
RIvNon |> filter(p.value < 0.01, estimate < 0) # okay!
RIvNon |> filter(p.value < 0.01, estimate > 0) # whew! probabilities positively correlated in *all* years

median(RIvNon$estimate)

ggplot(RIvNon, aes(x=estimate)) + geom_histogram() # and yeah the skew is in the direction we want


# YUBR --------------------------------
yubr.histStack <- raster::stack(sapply(yubr.files, function(x) crop(raster::raster(x), MojExt)))
names(yubr.histStack) <- paste("prFL",1900:2022,sep=".")
projection(yubr.histStack)<-CRS("+init=epsg:4269")

yubr.histStack

yubr.maskHist <- mask(yubr.histStack, st_transform(yubr.pres[,2], crs=4269))

yubr.hist.flowering <- cbind(coordinates(yubr.maskHist), as.data.frame(yubr.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(yubr.hist.flowering)

write.table(yubr.hist.flowering, "output/historic_flowering_reconst_YUBR.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


# YUJA --------------------------------
yuja.histStack <- raster::stack(sapply(yuja.files, function(x) crop(raster::raster(x), MojExt)))
names(yuja.histStack) <- paste("prFL",1900:2022,sep=".")
projection(yuja.histStack)<-CRS("+init=epsg:4269")

yuja.histStack

yuja.maskHist <- mask(yuja.histStack, st_transform(yuja.pres[,2], crs=4269))

yuja.hist.flowering <- cbind(coordinates(yuja.maskHist), as.data.frame(yuja.maskHist)) %>% filter(!is.na(prFL.2009)) %>% rename(lon=x, lat=y) %>% pivot_longer(starts_with("prFL"), names_to="year", values_to="prFL") %>% mutate(year=as.numeric(gsub("prFL\\.(\\d+)", "\\1", year)))

glimpse(yuja.hist.flowering)

write.table(yuja.hist.flowering, "output/historic_flowering_reconst_YUJA.csv", sep=",", col.names=TRUE, row.names=FALSE, quote=FALSE)


# YUJA+YUBR vs rangewide base model --------------------------------

# YUBR vs rangewide
yubr.mod <- read_rds("output/BART/bart.model.yubr.rds")
summary(yubr.mod) # best cutoff is 0.28

yubr.compare <- yubr.hist.flowering |> rename(yubr.prFL=prFL) |> mutate(lat=as.character(round(lat,3)), lon=as.character(round(lon,3)), year=as.character(year)) |> left_join(jotr.hist.flowering |> mutate(lat=as.character(round(lat,3)), lon=as.character(round(lon,3)), year=as.character(year))) |> mutate(lat=as.numeric(lat), lon=as.numeric(lon), year=as.numeric(year))
glimpse(yubr.compare) # okay that worked, ack

yubrtest <- yubr.compare |> filter(year>=2008)

yubrVbase <- table(yubr.compare$prFL>0.26, yubr.compare$yubr.prFL>0.28) # hmmmm!

sum(yubrVbase[c(1,4)])/sum(yubrVbase) # overall 0.62 agreed

chisq.test(yubrVbase) # p < 2.2e-16, natch

yubrVjotr <- yubrtest |> group_by(year) |> do(broom::tidy(cor.test(~prFL+yubr.prFL, data=., method="spearman"))) %>% ungroup()

glimpse(yubrVjotr)
yubrVjotr |> filter(p.value < 0.01, estimate < 0) # okay!
yubrVjotr |> filter(p.value < 0.01, estimate > 0) # whew! probabilities positively correlated in *all* years

min(yubrVjotr$estimate)
median(yubrVjotr$estimate)


# YUJA vs rangewide
yuja.mod <- read_rds("output/BART/bart.model.yuja.rds")
summary(yuja.mod) # best cutoff is 0.16

yuja.compare <- yuja.hist.flowering |> rename(yuja.prFL=prFL) |> mutate(lat=as.character(round(lat,3)), lon=as.character(round(lon,3)), year=as.character(year)) |> left_join(jotr.hist.flowering |> mutate(lat=as.character(round(lat,3)), lon=as.character(round(lon,3)), year=as.character(year))) |> mutate(lat=as.numeric(lat), lon=as.numeric(lon), year=as.numeric(year))
glimpse(yuja.compare) # okay that worked, ack

yujatest <- yuja.compare |> filter(year>=2008)

yujaVbase <- table(yuja.compare$prFL>0.26, yuja.compare$yuja.prFL>0.16) # hmmmm!

sum(yujaVbase[c(1,4)])/sum(yujaVbase) # overall 0.67 agreed

chisq.test(yujaVbase) # p < 2.2e-16, natch

yujaVjotr <- yujatest |> group_by(year) |> do(broom::tidy(cor.test(~prFL+yuja.prFL, data=., method="spearman"))) %>% ungroup()

glimpse(yujaVjotr)
yujaVjotr |> filter(p.value < 0.01, estimate < 0) # okay!
yujaVjotr |> filter(estimate < 0) # okay!
yujaVjotr |> filter(p.value < 0.01, estimate > 0) # whew! probabilities positively correlated in *all* years

min(yujaVjotr$estimate)
median(yujaVjotr$estimate)


#-------------------------------------------------------------------------
# Validation observations

# best-power prediction thresholds from the original models
# jotr: 0.26
# jotr RI: 0.25

# YUBR: 0.28
# YUJA: 0.16

glimpse(vobs) # 167 records

valid <- NULL # initialize, in a lazy way

for(yr in unique(vobs$year)){

# yr <- 2005

vsub <- filter(vobs, year==yr)

vsub$jotr_prFlr <- raster::extract(jotr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$jotr_Flr <- raster::extract(jotr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.26

vsub$jotr_RI.prFlr <- raster::extract(jotr.ri.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$jotr_RI.Flr <- raster::extract(jotr.ri.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.25


vsub$YUBR_prFlr <- raster::extract(yubr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$YUBR_Flr <- raster::extract(yubr.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.28

vsub$YUJA_prFlr <- raster::extract(yuja.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")])
vsub$YUJA_Flr <- raster::extract(yuja.histStack[[paste("prFL.",yr,sep="")]], vsub[,c("lon", "lat")]) >= 0.16

valid <- rbind(vsub, valid)

}

valid <- filter(valid, !is.na(jotr_prFlr))

glimpse(valid)
table(valid$obs_by)
head(valid) # cool
tail(valid)

write.table(valid, "output/expert_obs_validations.csv", sep=",", col.names=TRUE, row.names=FALSE)


# RANGE-WIDE -------------

# without RI
performance(prediction(valid$jotr_prFlr, valid$obs_flowers), "auc")@y.values[[1]] # 0.60 overall

t.test(jotr_prFlr~obs_flowers, data=valid) # p = 0.03 okay!

# jotr, with RI
performance(prediction(valid$jotr_RI.prFlr, valid$obs_flowers), "auc")@y.values[[1]] # 0.56 overall

t.test(jotr_RI.prFlr~obs_flowers, data=valid) # p = 0.11 hm

# SINGLE SPECIES -------------
vYUBR <- filter(valid, type=="YUBR")
vYUJA <- filter(valid, type=="YUJA")

# YUBR
performance(prediction(vYUBR$YUBR_prFlr, vYUBR$obs_flowers), "auc")@y.values[[1]] # 0.60 overall, whee
performance(prediction(valid$YUBR_prFlr, valid$obs_flowers), "auc")@y.values[[1]] # 0.51 overall, whee

t.test(YUBR_prFlr~obs_flowers, data=vYUBR) # p = 0.1

# YUJA
performance(prediction(vYUJA$YUJA_prFlr, vYUJA$obs_flowers), "auc")@y.values[[1]] # 0.48 overall, whee
performance(prediction(valid$YUJA_prFlr, valid$obs_flowers), "auc")@y.values[[1]] # 0.57 overall, whee

t.test(YUJA_prFlr~obs_flowers, data=vYUJA) # n.s. AND wrong direction?



# visualize accuracy with original data and validations ------------------

jotr.mod <- read_rds("output/BART/bart.model.Jotr.rds")
names(jotr.mod)

obs_valid <- summary(jotr.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data")
glimpse(obs_valid) 

exp_valid <- valid |> dplyr::select(jotr_prFlr, obs_flowers) |> rename(fitted=jotr_prFlr, observed=obs_flowers) |> mutate(observed=as.numeric(observed), type="Validation records")
glimpse(exp_valid)

accu <- rbind(obs_valid, exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.2635843)

# model accuracy figure
mod_accuracy <- ggplot(accu, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.26, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + facet_wrap("type", nrow=2) + 

annotate(geom="text", x=0.27, y=1.45, label="Classfication cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

# and then also the predictor selection figure
jotr.varimp <- read_rds("output/BART/bart.varimp.Jotr.rds")

jotr.varimp$data <- jotr.varimp$data |> mutate(trees = factor(trees, c(10,20,50,100,150,200)))

levels(jotr.varimp$data$names) <- c("Delta[Y1-2]*PPT", "Delta[Y0-1]*PPT", "Max*VPD[Y0]", "Delta[Y0-1]*Min*VPD", "Min*Temp[Y0]", "Delta[Y0-1]*Max*Temp", "PPT[Y0]", "PPT[Y1]", "Min*VPD[Y0]", "Delta[Y0-1]*Max*VPD", "Max*Temp[Y0]", "PPT[Y2]", "Delta[Y0-1]*Min*Temp")

jotr.varimp$labels$group <- "Trees"
jotr.varimp$labels$colour <- "Trees"

label_parse <- function(breaks){ parse(text=breaks) } # need this, for reasons

mod_varimp <- ggplot(jotr.varimp$data, aes(x=names, y=imp, color=trees, group=trees)) + geom_line(linewidth=0.75) + geom_point(size=1.5, color="gray30") + scale_color_manual(values=c('#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84'), name="Trees") + labs(y="Relative contribution") + scale_x_discrete(label=label_parse) + theme_minimal() + theme(legend.position=c(0.8, 0.75), axis.text=element_text(size=9), axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank(), legend.text=element_text(size=8), legend.title=element_text(size=8), legend.key.size=unit(0.15, "in"))  # okay nice


p <- read_rds("output/BART/bart.model.Jotr.partials.rds")

partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

partplot <- ggplot(partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))


{cairo_pdf("output/figures/Fig02_predictors_accuracy_partials.pdf", width=8, height=6)

ggdraw() + draw_plot(mod_varimp, 0, 0.53, 0.45, 0.44) + draw_plot(mod_accuracy, 0, 0, 0.45, 0.55) + draw_plot(partplot, 0.45, 0, 0.55, 1) + draw_plot_label(label=c("A", "B", "C"), x=c(0, 0, 0.45), y=c(1,0.55,1))

}
dev.off()


# figure to do the above but with the two species separately -------------

# first YUJA 
yuja.mod <- read_rds("output/BART/bart.model.yuja.rds")
names(yuja.mod)

yuja_obs_valid <- summary(yuja.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data")
glimpse(yuja_obs_valid) 

yuja_exp_valid <- valid |> filter(type=="YUJA") |> dplyr::select(YUJA_prFlr, obs_flowers) |> rename(fitted=YUJA_prFlr, observed=obs_flowers) |> mutate(observed=as.numeric(observed), type="Validation records")
glimpse(yuja_exp_valid)

yuja_accu <- rbind(yuja_obs_valid, yuja_exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.1589413)

# model accuracy figure
yuja_accuracy <- ggplot(yuja_accu, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.16, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + facet_wrap("type", nrow=2) + 

annotate(geom="text", x=0.17, y=1.45, label="Classfication cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

labs(title="Model accuracy, YUJA") + theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

yuja_accuracy

yuja_p <- partial(yuja.mod, jotr.preds, trace=FALSE, smooth=5)

yuja_partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", yuja_p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", yuja_p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", yuja_p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", yuja_p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", yuja_p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", yuja_p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

yuja_partplot <- ggplot(yuja_partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)", title="Partial effects, YUJA") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

yuja_partplot


# now YUBR
yubr.mod <- read_rds("output/BART/bart.model.yubr.rds")
names(yubr.mod)

yubr_obs_valid <- summary(yubr.mod)$data |> dplyr::select(fitted, observed) |> mutate(type="Training data")
glimpse(yubr_obs_valid) 

yubr_exp_valid <- valid |> filter(type=="YUBR") |> dplyr::select(YUBR_prFlr, obs_flowers) |> rename(fitted=YUBR_prFlr, observed=obs_flowers) |> mutate(observed=as.numeric(observed), type="Validation records")
glimpse(yubr_exp_valid)

yubr_accu <- rbind(yubr_obs_valid, yubr_exp_valid) |> mutate(type=factor(type, c("Training data", "Validation records")), classified=fitted>0.2776204)

# model accuracy figure
yubr_accuracy <- ggplot(yubr_accu, aes(y=as.factor(observed), x=fitted)) + geom_vline(xintercept=0.28, linewidth=0.5, color="gray60") + geom_jitter(height=0.15, aes(fill=classified, color=classified, shape=classified, size=classified, alpha=classified)) + geom_boxplot(fill=NA, color="black", linewidth=0.5, shape=20, outlier.color=NA, width=0.5) + facet_wrap("type", nrow=2) + 

annotate(geom="text", x=0.29, y=1.45, label="Classfication cutoff", hjust=0, color="gray40", size=3) +

scale_y_discrete(labels=c("False", "True")) + labs(y="Flowers observed", x="Predicted Pr(Flowers)") +

scale_fill_manual(values=c("#ffffcc", NA)) +
scale_color_manual(values=c("gray60", "#253494")) +
scale_shape_manual(values=c(21, 19)) + 
scale_size_manual(values=c(1, 1.2)) + 
scale_alpha_manual(values=c(0.75, 0.5)) + 

labs(title="Model accuracy, YUBR") + theme_minimal(base_size=10) + theme(legend.position="none", plot.margin=unit(c(0.1,0.1,0.1,0.1),"in"), strip.text=element_text(size=10, hjust=0))

yubr_accuracy

yubr_p <- partial(yubr.mod, jotr.preds, trace=FALSE, smooth=5)

yubr_partvals <- rbind(
				data.frame(predictor="Delta[Y1-2]*PPT~(mm)", yubr_p[[1]]$data),
				data.frame(predictor="Delta[Y0-1]*PPT~(mm)", yubr_p[[2]]$data),
				data.frame(predictor="Max*VPD[Y0]~(hPa)", yubr_p[[3]]$data),
				data.frame(predictor="Delta[Y0-1]*Min*VPD~(hPa)", yubr_p[[4]]$data),
				data.frame(predictor="Min*Temp[Y0]~(degree*C)", yubr_p[[5]]$data),
				data.frame(predictor="Delta[Y0-1]*Max*Temp~(degree*C)", yubr_p[[6]]$data)
				) |> mutate(predictor=factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

yubr_partplot <- ggplot(yubr_partvals) + geom_ribbon(aes(x=x, ymin=q05, ymax=q95), fill="#41b6c4") + geom_line(aes(x=x, y=med), color="white") + facet_wrap("predictor", nrow=3, labeller="label_parsed", scale="free") + labs(y="Marginal Pr(Flowers)", title="Partial effects, YUBR") + theme_minimal() + theme(axis.title.x=element_blank(), panel.spacing=unit(0.2,"in"))

yubr_partplot


{cairo_pdf("output/figures/SIFig_subspecies_accuracy_partials.pdf", width=6.5, height=8)

ggdraw() + draw_plot(yubr_accuracy, 0, 0.6, 0.5, 0.4) + draw_plot(yuja_accuracy, 0.5, 0.6, 0.5, 0.4) + draw_plot(yubr_partplot, 0, 0, 0.5, 0.6) + draw_plot(yuja_partplot, 0.5, 0, 0.5, 0.6) + draw_plot_label(label=c("A", "B", "C", "D"), x=c(0, 0, 0.5, 0.5), y=c(1, 0.6, 1, 0.6))

}
dev.off()




#-------------------------------------------------------------------------
# TRENDS within cells

# cor(flr, yr) --- crude but it's a start -----------------

# not sure this is the best approach?
cellcors <- jotr.hist.flowering %>% group_by(lon,lat) %>% do(broom::tidy(cor.test(~prFL+year, data=., method="spearman"))) %>% ungroup() # rangewide model

glimpse(cellcors)

quantile(cellcors$estimate, c(0.05, 0.5, 0.975))

{cairo_pdf("output/figures/prFL-vs-time_correlations_jotr.pdf", width=3, height=2.5)

ggplot(cellcors, aes(x=estimate, fill=p.value<=0.01, color=p.value<=0.01)) + geom_histogram(color=NA) + 

scale_fill_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +
scale_color_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +

#annotate("text", x=0.1, y=150, label="All cells", color=park_palette("JoshuaTree")[7], fontface="bold", family="Arial Narrow") + 

#annotate("text", x=-0.15, y=50, label="Cells with\ncorrelation\np < 0.01", color=park_palette("JoshuaTree")[2], fontface="bold", family="Arial Narrow", lineheight=0.8) + 

geom_vline(xintercept=0) +

labs(x=expression("Spearman's"~rho), y="Grid cells") + theme_minimal() + theme(legend.position="none", plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()

cellcorsRI <- jotr.hist.flowering %>% group_by(lon,lat) %>% do(broom::tidy(cor.test(~RI.prFL+year, data=., method="spearman"))) %>% ungroup() # rangewide model

glimpse(cellcorsRI)

{cairo_pdf("output/figures/RI.prFL-vs-time_correlations_jotr.pdf", width=3, height=2.5)

ggplot(cellcorsRI, aes(x=estimate, fill=p.value<=0.01, color=p.value<=0.01)) + geom_histogram(color=NA) + 

scale_fill_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +
scale_color_manual(values=park_palette("JoshuaTree")[c(7, 2)]) +

#annotate("text", x=0.1, y=150, label="All cells", color=park_palette("JoshuaTree")[7], fontface="bold", family="Arial Narrow") + 

#annotate("text", x=-0.15, y=50, label="Cells with\ncorrelation\np < 0.01", color=park_palette("JoshuaTree")[2], fontface="bold", family="Arial Narrow", lineheight=0.8) + 

geom_vline(xintercept=0) +

labs(x=expression("Spearman's"~rho), y="Grid cells") + theme_minimal() + theme(legend.position="none", plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()


#-------------------------------------------------------------------------
# Flowering years

# best-power prediction thresholds from the original models
# jotr: 0.26
# jotr with RI: 0.25

# YUBR: 0.30
# YUJA: 0.23

flyrs.jotr <- jotr.hist.flowering %>% group_by(lon,lat) %>% 
summarize(
	flyrs_all=length(which(prFL>=0.26)), 
	flyrs_1900_1929=length(which(prFL>=0.26 & year<=1929)), 
	flyrs_1990_2019=length(which(prFL>=0.26 & year>=1990 & year<=2019)),
	flyrs_RI_all=length(which(RI.prFL>=0.25)), 
	flyrs_RI_1900_1929=length(which(RI.prFL>=0.25 & year<=1929)), 
	flyrs_RI_1990_2019=length(which(RI.prFL>=0.25 & year>=1990 & year<=2019)),
	) |> mutate(
	flyrs_change=flyrs_1990_2019-flyrs_1900_1929,
	flyrs_RI_change=flyrs_RI_1990_2019-flyrs_RI_1900_1929
	) # rangewide model
range(jotr.hist.flowering$year) # remember: 1900-2022
glimpse(flyrs.jotr)

write.table(flyrs.jotr, "output/jotr_reconstructed_flowering_years.csv", sep=",", col.names=TRUE, row.names=FALSE)
# flyrs.jotr <- read.csv("output/jotr_reconstructed_flowering_years.csv")

# without RI
quantile(flyrs.jotr$flyrs_all, c(0.025,0.5,0.975))
quantile(flyrs.jotr$flyrs_all, c(0.025,0.5,0.975))/123 # median 0.30

quantile(flyrs.jotr$flyrs_1900_1929, c(0.025,0.5,0.975))/30 # median 0.27
quantile(flyrs.jotr$flyrs_1990_2019, c(0.025,0.5,0.975))/30 # median 0.33

quantile(flyrs.jotr$flyrs_change, c(0.025,0.5,0.975)) # median 2


# with RI
quantile(flyrs.jotr$flyrs_RI_all, c(0.025,0.5,0.975))
quantile(flyrs.jotr$flyrs_RI_all, c(0.025,0.5,0.975))/123 # median 0.24 (!)

quantile(flyrs.jotr$flyrs_RI_1900_1929, c(0.025,0.5,0.975))/30 # median 0.23
quantile(flyrs.jotr$flyrs_RI_1990_2019, c(0.025,0.5,0.975))/30 #  median 0.23 ohhh fascinating

quantile(flyrs.jotr$flyrs_RI_change, c(0.025,0.5,0.975)) # median 1


ggplot(flyrs.jotr, aes(x=flyrs_1900_1929)) + geom_histogram(bins=30)
ggplot(flyrs.jotr, aes(x=flyrs_1990_2019)) + geom_histogram(bins=30)


{cairo_pdf("output/figures/flowering-years_jotr.pdf", width=3.5, height=2.5)

ggplot(flyrs.jotr, aes(x=flyrs_all)) + geom_histogram(fill="#ccece6") + labs(x="Projected flowering years, 1900-2022", y="Grid cells") + 

geom_vline(xintercept=median(flyrs.jotr$flyrs_all), color="#006d2c") +

theme_minimal() + theme(legend.position="none", plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()

{cairo_pdf("output/figures/RI_flowering-years_jotr.pdf", width=3.5, height=2.5)

ggplot(flyrs.jotr, aes(x=flyrs_RI_all)) + geom_histogram(fill="#ccece6") + labs(x="Projected flowering years, 1900-2022", y="Grid cells") + 

geom_vline(xintercept=median(flyrs.jotr$flyrs_RI_all), color="#006d2c") +

theme_minimal() + theme(legend.position="none", plot.margin=margin(0.1,0.15,0.1,0.15, "inches"))

}
dev.off()


# Flowering years, mapped
{cairo_pdf("output/figures/flyrs_map_jotr.pdf", width=6, height=5)

ggplot(flyrs.jotr, aes(x=lon, y=lat, fill=flyrs_all)) + geom_tile() + 

#coord_fixed() + 

scale_fill_distiller(type="seq", palette=2, direction=1, name="Years flowering predicted") + labs(x="Longitude", y="Latitude", title="Predicted flowering years, 1900-2022") + 

theme_minimal() + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.15, "inches"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.2,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6))

}
dev.off()

#-------------------------------------------------------------------------
# Change in flowering years, trends in key predictors

ggplot(flyrs.jotr, aes(x=flyrs_1990_2019-flyrs_1900_1929)) + geom_histogram(bins=20)
ggplot(flyrs.jotr, aes(x=flyrs_RI_1990_2019-flyrs_RI_1900_1929)) + geom_histogram(bins=20)

DelFlyrs <- flyrs.jotr |> mutate(FlyrsBase = flyrs_1990_2019-flyrs_1900_1929, FlyrsRI = flyrs_RI_1990_2019-flyrs_RI_1900_1929) |> dplyr::select(lon, lat, FlyrsBase, FlyrsRI) |> pivot_longer(contains("Flyrs"), names_to="model", values_to="FlyrsChange") |> mutate(model=factor(model, c("FlyrsBase", "FlyrsRI")))

levels(DelFlyrs$model) <- c("Base Model", "RI Model")


# map context
states <- read_sf(dsn = "../data/spatial/10m_cultural/", lay= "ne_10m_admin_1_states_provinces")
coast <- read_sf("../data/spatial/10m_physical/ne_10m_coastline", "ne_10m_coastline")


# change in flowering years, baseline and RI model
flyrsChange <- ggplot() + 
	geom_sf(data=coast, color="slategray2", linewidth=2.5) + 
	geom_sf(data=states, fill="cornsilk3", color="antiquewhite4") + 
	
	geom_tile(data=DelFlyrs, aes(x=lon, y=lat, fill=FlyrsChange)) + 
	
	scale_fill_distiller(type="div", palette="PRGn", direction=1, name="Change in flowering years,\n1990-2019 vs 1900-1929") + labs(x="Longitude", y="Latitude") + 
	
	facet_wrap("model", nrow=2) +
	
	coord_sf(xlim = c(-119, -112.5), ylim = c(33.55, 37.75), expand = FALSE) +
	theme_minimal(base_size=8) + theme(legend.position="bottom", legend.key.width=unit(0.15, "inches"), legend.key.height=unit(0.1, "in"), axis.text=element_blank(), axis.title=element_blank(), plot.margin=unit(c(0.1,0.1,0.1,0.1), "inches"), panel.spacing=unit(0.1,"inches"), legend.spacing.y=unit(0.1,"inches"), legend.box="horizontal", legend.text=element_text(size=6), legend.title=element_text(size=8), strip.text=element_text(size=10), panel.background=element_rect(fill="slategray3", color="black"), panel.grid=element_blank())

flyrsChange


# change in predictors .................................
sample_sites <- flyrs.jotr |> dplyr::select(lat, lon)

sample_predictor_history <- data.frame(matrix(0,0,3+length(jotr.preds)))
names(sample_predictor_history) <- c("lat", "lon", "year", jotr.preds)

for(yr in 1900:2022){ # loop over years

sample_predictor_history <- rbind(sample_predictor_history, data.frame(sample_sites,year=yr, raster::extract(brick(paste("data/PRISM/derived_predictors/PRISM_derived_predictors_",yr,".gri",sep="")), sample_sites[,c("lon","lat")])[,jotr.preds]))

}

# now compile trends ...
pred_history_long <- sample_predictor_history |> pivot_longer(all_of(jotr.preds), names_to="predictor", values_to="value") |> mutate(predictor = factor(predictor, jotr.preds))

levels(pred_history_long$predictor) <- c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")


# change in MEAN from 1900-1929 to 1990-2019
predictor_mean_early <- pred_history_long |> filter(year%in%1900:1929) |> group_by(lat, lon, predictor) |> summarize(mean_1900_1929=mean(value))

predictor_mean_late <- pred_history_long |> filter(year%in%1990:2019) |> group_by(lat, lon, predictor) |> summarize(mean_1990_2019=mean(value))

predictor_mean_change <- full_join(predictor_mean_early, predictor_mean_late) |> mutate(mean_change=mean_1990_2019-mean_1900_1929) |> pivot_longer(starts_with("mean"), values_to="value", names_to="mean_stat")
glimpse(predictor_mean_change)

# BART predictor change in means
pred_change <- ggplot(filter(predictor_mean_change, mean_stat!="mean_change"), aes(x=value, fill=mean_stat)) + geom_histogram(position="dodge") + facet_wrap("predictor", scale="free", labeller="label_parsed", nrow=3) + 

scale_fill_manual(values=park_palette("JoshuaTree")[c(1,2)], labels=c("1900-1929", "1990-2019"), name="Mean for") + 

labs(y="Grid cells", x="Mean value") + theme_minimal(base_size=9) + theme(legend.position="bottom", legend.key.size=unit(0.1, "in")) 



# change in flowering years vs other factors ...
glimpse(flyrs.jotr)

jotr.flr.change <- cbind(flyrs.jotr, raster::extract(raster("../data/Yucca/ybrev_elev.grd"), flyrs.jotr[,c("lon","lat")])) 
colnames(jotr.flr.change)[11] <- "elev_m"

glimpse(jotr.flr.change)

cor.test(~flyrs_change+lat, data=jotr.flr.change, method="sp")
cor.test(~flyrs_RI_change+lat, data=jotr.flr.change, method="sp") # huh

cor.test(~flyrs_change+elev_m, data=jotr.flr.change, method="sp")
cor.test(~flyrs_RI_change+elev_m, data=jotr.flr.change, method="sp") # okay, largely aligned here


# paired predictor and flowering year change ...............
predictor_change <- predictor_mean_change |> mutate(mean_stat = factor(mean_stat, levels=c("mean_1900_1929", "mean_1990_2019", "mean_change"))) |> rename(timeframe=mean_stat, pred_value=value)
levels(predictor_change$timeframe) <- c("from_1900_1929", "from_1990_2019", "change")

flyrs_change <- flyrs.jotr |> pivot_longer(starts_with("flyrs_"), names_to="timeframe", values_to="flyrs") |> mutate(ri.model=grepl("RI", timeframe)) 

flyrs_change$timeframe[grepl("change", flyrs_change$timeframe)] <- "change"
flyrs_change$timeframe[grepl("1900_1929", flyrs_change$timeframe)] <- "from_1900_1929"
flyrs_change$timeframe[grepl("1990_2019", flyrs_change$timeframe)] <- "from_1990_2019"
flyrs_change$timeframe[grepl("all", flyrs_change$timeframe)] <- "from_all"
table(flyrs_change$timeframe)
table(predictor_change$timeframe)

flyrs_pred_change <- flyrs_change |> filter(timeframe!="from_all") |> mutate(timeframe=factor(timeframe, c("from_1900_1929", "from_1990_2019", "change"))) |> left_join(predictor_change)

glimpse(flyrs_pred_change)
table(flyrs_pred_change$timeframe)

write.table(flyrs_pred_change, "output/jotr_flowering_predictors_change.csv", sep=",", col.names=TRUE, row.names=FALSE)

# flyrs_pred_change <- read.csv("output/jotr_flowering_predictors_change.csv") |> mutate(predictor = factor(predictor, c("Delta[Y1-2]*PPT~(mm)", "Delta[Y0-1]*PPT~(mm)", "Max*VPD[Y0]~(hPa)", "Delta[Y0-1]*Min*VPD~(hPa)", "Min*Temp[Y0]~(degree*C)", "Delta[Y0-1]*Max*Temp~(degree*C)")))

ggplot(jotr.flr.change, aes(x=elev_m, y=flyrs_change)) + geom_point(alpha=0.2) + geom_smooth(method="lm")


# BART predictor change vs flowering changes

plot_slice <- filter(flyrs_pred_change, timeframe!="change", !ri.model) |> group_by(timeframe, predictor, ri.model) |> slice_sample(n=300)

plot_slice_summ <- plot_slice |> group_by(timeframe, predictor) |> summarize(med_pred = median(pred_value), med_flyrs = median(flyrs)) |> pivot_wider(names_from=timeframe, values_from=c(med_pred, med_flyrs))


pred_flr_change <- ggplot(plot_slice, aes(x=pred_value, y=flyrs, color=timeframe)) + geom_point(alpha=0.25, size=1.25) + geom_segment(data=plot_slice_summ, aes(x = med_pred_from_1900_1929, y = med_flyrs_from_1900_1929, xend = med_pred_from_1990_2019, yend = med_flyrs_from_1990_2019), arrow = arrow(length = unit(0.15, "cm")), color="black") +

facet_wrap("predictor", scale="free", labeller="label_parsed", nrow=3) + 

scale_color_manual(values=c("#0571b0","#e66101"), labels=c("1900-1929", "1990-2019"), name="Value for") + 

labs(y="Flowering years predicted by base model", x="Mean predictor value") + theme_minimal(base_size=9) + theme(legend.position="bottom", legend.key.size=unit(0.1, "in"))


# write out composite figure

{cairo_pdf(file="output/figures/Fig04_flyrs_predictors_change.pdf", width=6, height=5)

ggdraw() + draw_plot(flyrsChange, 0, 0, 0.44, 1) + draw_plot(pred_flr_change, 0.43, 0, 0.53, 1) + draw_plot_label(label=c("A", "B"), x=c(0,0.43), y=1)

}
dev.off()
