rm(list = ls())
library(plyr);library(dplyr);library(ggplot2)

## budburst groups: budburst & no budburst; single tree base

################################################################################
## Planetscope
load("out/planetscope/ndvi_long_format_phenoclasses_planetscope.RData")
planetscope <- ndvi_long_pheno_planetscope; rm(ndvi_long_pheno_planetscope)

# checking for outliers within budburst groups per tree
boxplot(planetscope$ndvi[which(planetscope$budburst == T)]) #seemingly no outliers
boxplot(planetscope$ndvi[which(planetscope$budburst == F)]) #outliers

# remove outliers for "no budburst" based the IQR 
planetscope <- planetscope[-which(planetscope$budburst == F & planetscope$ndvi %in% c(boxplot.stats(planetscope$ndvi[which(planetscope$budburst == F)])$out)),]

################################################################################
## Sentinel-2
load("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.RData")
sentinel2 <- ndvi_long_pheno_sentinel2; rm(ndvi_long_pheno_sentinel2)
sentinel2$date <- as.Date(sentinel2$date)

# checking for outliers within budburst groups per tree
boxplot(sentinel2$ndvi[which(sentinel2$budburst == T)]) #one outliers
boxplot(sentinel2$ndvi[which(sentinel2$budburst == F)]) #several outliers

#remove outliers for "budburst" based the IQR 
sentinel2 <- sentinel2[-which(sentinel2$budburst == T & sentinel2$ndvi %in% c(boxplot.stats(sentinel2$ndvi[which(sentinel2$budburst == T)])$out)),]

#remove outliers for "no budburst" based the IQR 
sentinel2 <- sentinel2[-which(sentinel2$budburst == F & sentinel2$ndvi %in% c(boxplot.stats(sentinel2$ndvi[which(sentinel2$budburst == F)])$out)),]


################################################################################
## TreeTalker
load("out/treetalker/ndvi_long_format_phenoclasses_treetalker.RData")
treetalker <- ndvi_long_format_phenoclasses_treetalker; rm(ndvi_long_format_phenoclasses_treetalker)

# checking for outliers within budburst groups per tree
boxplot(treetalker$ndvi[which(treetalker$budburst == T)]) #sseveral outliers
boxplot(treetalker$ndvi[which(treetalker$budburst == F)]) #sseveral outliers

#remove outliers for "budburst" based the IQR 
treetalker <- treetalker[-which(treetalker$budburst == T & treetalker$ndvi %in% c(boxplot.stats(treetalker$ndvi[which(treetalker$budburst == T)])$out)),]

#remove outliers for "no budburst" based the IQR 
treetalker <- treetalker[-which(treetalker$budburst == F & treetalker$ndvi %in% c(boxplot.stats(treetalker$ndvi[which(treetalker$budburst == F)])$out)),]
