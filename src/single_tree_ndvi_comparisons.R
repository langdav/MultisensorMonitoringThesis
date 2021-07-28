rm(list = ls())
library(plyr);library(dplyr);library(ggplot2)

# compare daily per tree NDVI on daily base

################################################################################
## Planetscope
load("out/planetscope/ndvi_long_format_phenoclasses_planetscope.RData")
planetscope <- ndvi_long_pheno_planetscope; rm(ndvi_long_pheno_planetscope)

# checking for outliers on per tree and per day base before calculating daily means
planetscope_mean <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(planetscope_mean) <- c(colnames(planetscope)[1:2],"ndvi_mean", "ndvi_sd", colnames(planetscope)[4:5])
for(tree in unique(planetscope$tree_id)){
  days <- unique(planetscope$date[which(planetscope$tree_id == tree)])
  for(day in 1:length(days)){
    ndvi <- planetscope$ndvi[which(planetscope$tree_id == tree & planetscope$date == days[day])]
    outliers <- which(ndvi %in% c(boxplot.stats(ndvi)$out)) #locate outliers
    
    if(length(outliers) != 0) {
      ndvi_mean <- mean(ndvi[-outliers], na.rm = T)
      ndvi_sd <- sd(ndvi[-outliers], na.rm = T)
    } else {
      ndvi_mean <- mean(ndvi, na.rm = T)
      ndvi_sd <- sd(ndvi, na.rm = T)
    }
    
    bb <- unique(planetscope$budburst[which(planetscope$tree_id == tree & planetscope$date == days[day])])
    bbpc <- unique(planetscope$budburst_perc[which(planetscope$tree_id == tree & planetscope$date == days[day])])
    
    df_tmp <- data.frame(tree, days[day], ndvi_mean, ndvi_sd, bb, bbpc)
    colnames(df_tmp) <- colnames(planetscope_mean)
    
    planetscope_mean <- rbind(planetscope_mean, df_tmp)
  }
}


################################################################################
## Sentinel-2
load("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.RData")
sentinel2 <- ndvi_long_pheno_sentinel2; rm(ndvi_long_pheno_sentinel2)
sentinel2$date <- as.Date(sentinel2$date)


################################################################################
## TreeTalker
load("out/treetalker/ndvi_long_format_phenoclasses_treetalker.RData")
treetalker <- ndvi_long_format_phenoclasses_treetalker; rm(ndvi_long_format_phenoclasses_treetalker)

# checking for outliers on per tree and per day base before calculating daily means
treetalker_mean <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(treetalker_mean) <- c(colnames(treetalker)[1:2],"ndvi_mean", "ndvi_sd", colnames(treetalker)[5:6])
for(tree in unique(treetalker$tree_id)){
  days <- unique(treetalker$date[which(treetalker$tree_id == tree)])
  for(day in 1:length(days)){
    ndvi <- treetalker$ndvi[which(treetalker$tree_id == tree & treetalker$date == days[day])]
    outliers <- which(ndvi %in% c(boxplot.stats(ndvi)$out)) #locate outliers
    
    if(length(outliers) != 0) {
      ndvi_mean <- mean(ndvi[-outliers], na.rm = T)
      ndvi_sd <- sd(ndvi[-outliers], na.rm = T)
    } else {
      ndvi_mean <- mean(ndvi, na.rm = T)
      ndvi_sd <- sd(ndvi, na.rm = T)
    }
    
    bb <- unique(treetalker$budburst[which(treetalker$tree_id == tree & treetalker$date == days[day])])
    bbpc <- unique(treetalker$budburst_perc[which(treetalker$tree_id == tree & treetalker$date == days[day])])
    
    df_tmp <- data.frame(tree, days[day], ndvi_mean, ndvi_sd, bb, bbpc)
    colnames(df_tmp) <- colnames(treetalker_mean)
    
    treetalker_mean <- rbind(treetalker_mean, df_tmp)
  }
}
