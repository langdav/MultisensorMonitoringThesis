#author: David Langenohl
#last modified: 25.08.2021
#description: create different all-in-one data frames
#NOTE:

rm(list = ls())
library(plyr);library(dplyr);library(rstatix);library(ggpubr)

# daily means per tree all-in-one
load("out/planetscope/outlier_free_daily_ndvi_mean_per_tree_with_phenoclasses_planetscope.RData");planetscope_mean <- as.data.frame(ungroup(planetscope_ndvi_mean_per_tree));rm(planetscope_ndvi_mean_per_tree)
load("out/sentinel2/ndvi_mean_per_tree_with_phenoclasses_sentinel2.RData");sen2_mean <- as.data.frame(ungroup(sen2_ndvi_mean_per_tree));rm(sen2_ndvi_mean_per_tree)
load("out/treetalker/outlier_free_ndvi_mean_per_tree_with_phenoclasses_treetalker.RData");treetalker_mean <- as.data.frame(ungroup(treetalker_ndvi_mean_per_tree));rm(treetalker_ndvi_mean_per_tree)
load("out/orthomosaic/outlier_free_ndvi_mean_per_tree_with_phenoclasses_orthomosaic.RData");orthomosaic_mean <- as.data.frame(ungroup(orthomosaic_ndvi_mean_per_tree));rm(orthomosaic_ndvi_mean_per_tree)

planetscope_mean$platform <- rep("planetscope", nrow(planetscope_mean))
sen2_mean$platform <- rep("sentinel2", nrow(sen2_mean))
treetalker_mean$platform <- rep("treetalker", nrow(treetalker_mean))
orthomosaic_mean$platform <- rep("orthomosaic", nrow(orthomosaic_mean))

aio_daily_ndvi_means <- rbind(planetscope_mean, sen2_mean, treetalker_mean, orthomosaic_mean)

# daily means over all trees per platform
aio_all_tree_means <- NULL
for(platform in unique(aio_daily_ndvi_means$platform)){
  cat("Processing", platform, "\n")
  days <- unique(aio_daily_ndvi_means$date[which(aio_daily_ndvi_means$platform == platform)])
  for(day in 1:length(days)){
    tmp <- aio_daily_ndvi_means %>% filter(platform == platform) %>% filter(date == days[day])
    m <- mean(tmp$ndvi_mean, na.rm = T)
    std <- sd(tmp$ndvi_mean, na.rm = T)
    bb <- ifelse(any(tmp$budburst)==T, T, F)
    
    aio_all_tree_means <- rbind(aio_all_tree_means, data.frame(platform = platform,
                                                               date = days[day],
                                                               doy = yday(days[day]),
                                                               ndvi_mean = m,
                                                               ndvi_sd = std,
                                                               budburst = bb))
  }
}

save(aio_daily_ndvi_means, file = "out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
save(aio_all_tree_means, file = "out/all_in_one/aio_daily_ndvi_all_trees_means.RData")

# daily medians per tree all-in-one
load("out/planetscope/outlier_free_daily_ndvi_median_per_tree_with_phenoclasses_planetscope.RData");planetscope_median <- as.data.frame(ungroup(planetscope_ndvi_median_per_tree));rm(planetscope_ndvi_median_per_tree)
load("out/sentinel2/ndvi_median_per_tree_with_phenoclasses_sentinel2.RData");sen2_median <- as.data.frame(ungroup(sen2_ndvi_median_per_tree));rm(sen2_ndvi_median_per_tree)
load("out/treetalker/outlier_free_ndvi_median_per_tree_with_phenoclasses_treetalker.RData");treetalker_median <- as.data.frame(ungroup(treetalker_ndvi_median_per_tree));rm(treetalker_ndvi_median_per_tree)
load("out/orthomosaic/outlier_free_ndvi_median_per_tree_with_phenoclasses_orthomosaic.RData");orthomosaic_median <- as.data.frame(ungroup(orthomosaic_ndvi_median_per_tree));rm(orthomosaic_ndvi_median_per_tree)

planetscope_median$platform <- rep("planetscope", nrow(planetscope_median))
sen2_median$platform <- rep("sentinel2", nrow(sen2_median))
treetalker_median$platform <- rep("treetalker", nrow(treetalker_median))
orthomosaic_median$platform <- rep("orthomosaic", nrow(orthomosaic_median))

aio_daily_ndvi_medians <- rbind(planetscope_median, sen2_median, treetalker_median, orthomosaic_median)

# daily means over all trees per platform
aio_all_tree_medians <- NULL
for(platform in unique(aio_daily_ndvi_medians$platform)){
  cat("Processing", platform, "\n")
  days <- unique(aio_daily_ndvi_medians$date[which(aio_daily_ndvi_medians$platform == platform)])
  for(day in 1:length(days)){
    tmp <- aio_daily_ndvi_medians %>% filter(platform == platform) %>% filter(date == days[day])
    m <- mean(tmp$ndvi_median, na.rm = T)
    std <- sd(tmp$ndvi_median, na.rm = T)
    bb <- ifelse(any(tmp$budburst)==T, T, F)
    
    aio_all_tree_medians <- rbind(aio_all_tree_medians, data.frame(platform = platform,
                                                               date = days[day],
                                                               doy = yday(days[day]),
                                                               ndvi_mean = m,
                                                               ndvi_sd = std,
                                                               budburst = bb))
  }
}

save(aio_daily_ndvi_medians, file = "out/all_in_one/aio_daily_ndvi_per_tree_medians.RData")
save(aio_all_tree_medians, file = "out/all_in_one/aio_daily_ndvi_all_trees_medians.RData")

# daily mean backscatter (Sentinel-1) over all trees
load("out/sentinel1/backscatter_all_with_phenoclasses_sentinel1.RData")
sen1_all_tree_means <- NULL
days <- unique(sen1_backscatter$date)
for(day in 1:length(days)){
  tmp <- sen1_backscatter %>% filter(date == days[day])
  m <- mean(tmp$sigma_ratio, na.rm = T)
  std <- sd(tmp$sigma_ratio, na.rm = T)
  bb <- ifelse(any(tmp$budburst)==T, T, F)
  
  sen1_all_tree_means <- rbind(sen1_all_tree_means, data.frame(platform = "sentinel1",
                                                             date = days[day],
                                                             doy = yday(days[day]),
                                                             sigma_ratio_mean = m,
                                                             sigma_ratio_sd = std,
                                                             budburst = bb))
}

save(sen1_all_tree_means, file = "out/all_in_one/sentinel1_daily_sigma_ratio_all_trees.RData")
