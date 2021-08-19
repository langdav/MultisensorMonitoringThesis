#author: David Langenohl
#last modified: 19.08.2021
#description: create different all-in-one data frames
#NOTE:

rm(list = ls())
library(plyr);library(dplyr);library(rstatix);library(ggpubr)

#load data
load("out/planetscope/outlier_free_daily_ndvi_mean_per_tree_with_phenoclasses_planetscope.RData");planetscope_mean <- as.data.frame(ungroup(planetscope_ndvi_mean_per_tree));rm(planetscope_ndvi_mean_per_tree)
load("out/sentinel2/ndvi_mean_per_tree_with_phenoclasses_sentinel2.RData");sen2_mean <- as.data.frame(ungroup(sen2_ndvi_mean_per_tree));rm(sen2_ndvi_mean_per_tree)
load("out/treetalker/outlier_free_ndvi_mean_per_tree_with_phenoclasses_treetalker.RData");treetalker_mean <- as.data.frame(ungroup(treetalker_ndvi_mean_per_tree));rm(treetalker_ndvi_mean_per_tree)
load("out/orthomosaic/outlier_free_ndvi_mean_per_tree_with_phenoclasses_orthomosaic.RData");orthomosaic_mean <- as.data.frame(ungroup(orthomosaic_ndvi_mean_per_tree));rm(orthomosaic_ndvi_mean_per_tree)

planetscope_mean$platform <- rep("planetscope", nrow(planetscope_mean))
sen2_mean$platform <- rep("sentinel2", nrow(sen2_mean))
treetalker_mean$platform <- rep("treetalker", nrow(treetalker_mean))
orthomosaic_mean$platform <- rep("orthomosaic", nrow(orthomosaic_mean))

aio_daily_ndvi_means <- rbind(planetscope_mean, sen2_mean, treetalker_mean, orthomosaic_mean)

# save resulting data frame
save(aio_daily_ndvi_means, file = "out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
# =============================================================

# --------------------------------------------------------------------------------------
# create NDVI means over all trees for each date ---------------------------------------

# load data
load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)

# detect and remove outliers
for(platform in unique(aio$platform)){
  days <- unique(aio$date[which(aio$platform == platform)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(aio$ndvi_mean[which(aio$platform == platform & 
                                                aio$date == days[day])])$out) != 0){
      aio <- aio[-which(aio$platform == platform & 
                          aio$date == days[day] & 
                          aio$ndvi_mean %in% c(boxplot.stats(aio$ndvi_mean[which(aio$platform == platform & 
                                                                                   aio$date == days[day])])$out)),]
    }
  }
}

# calculate means over all trees
aio_all_tree_means <- NULL
for(platform in unique(aio$platform)){
  cat("Processing", platform, "\n")
  days <- unique(aio$date[which(aio$platform == platform)])
  for(day in 1:length(days)){
    tmp <- aio[-which(aio$platform == platform & 
                        aio$date == days[day]),]
    m <- mean(tmp$ndvi_mean, na.rm = T)
    std <- sd(tmp$ndvi_mean, na.rm = T)
    bb <- ifelse(any(tmp$budburst)==T, T, F)
    
    aio_all_tree_means <- rbind(aio_all_tree_means, data.frame(platform = platform,
                                                               date = days[day],
                                                               ndvi_mean = m,
                                                               ndvi_sd = std,
                                                               budburst = bb))
  }
}

# save resulting data frame
save(aio_all_tree_means, file = "out/all_in_one/aio_all_tree_means_per_platform_and_day.RData")

