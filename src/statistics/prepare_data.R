rm(list = ls())
library(plyr);library(dplyr);library(rstatix);library(ggpubr)

################
## load data ##
##############
# Planetscope
load("out/planetscope/ndvi_long_format_phenoclasses_planetscope.RData")
planetscope <- ndvi_long_pheno_planetscope; rm(ndvi_long_pheno_planetscope)

# Sentinel-2
load("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.RData")
sentinel2 <- ndvi_long_pheno_sentinel2; rm(ndvi_long_pheno_sentinel2)
sentinel2$date <- as.Date(sentinel2$date)
sentinel2 <- rename(sentinel2, ndvi_mean = ndvi)

# Orthomosaic
load("out/orthomosaic/ndvi_long_format_phenoclasses_orthomosaic.RData") #only load if necessary, as resulting dataframe is huge and contains 8553222 rows
orthomosaic <- ndvi_long_pheno_orthomosaic; rm(ndvi_long_pheno_orthomosaic)

# TreeTalker
load("out/treetalker/ndvi_long_format_phenoclasses_treetalker.RData")
treetalker <- ndvi_long_format_phenoclasses_treetalker; rm(ndvi_long_format_phenoclasses_treetalker)

##########################################################################################
## calculate per tree daily mean NDVI values, if platform has multiple values per tree ##
########################################################################################
## Planetscope
# remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(planetscope$tree_id)){
  days <- unique(planetscope$date[which(planetscope$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(planetscope$ndvi[which(planetscope$tree_id == tree & 
                                                   planetscope$date == days[day])])$out) != 0){
      planetscope <- planetscope[-which(planetscope$tree_id == tree & 
                                          planetscope$date == days[day] & 
                                          planetscope$ndvi %in% c(boxplot.stats(planetscope$ndvi[which(planetscope$tree_id == tree & 
                                                                                                         planetscope$date == days[day])])$out)),]
    }
  }
}

# after removing outliers, calculate daily mean NDVI values
planetscope_mean_per_tree <- planetscope %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

# save resulting data frame
save(planetscope_mean_per_tree, file = "out/planetscope/outlier_free_daily_ndvi_mean_per_tree_planetscope.RData")

## Orthomosaic
# remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(orthomosaic$tree_id)){
  days <- unique(orthomosaic$date[which(orthomosaic$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(orthomosaic$ndvi[which(orthomosaic$tree_id == tree &
                                                   orthomosaic$date == days[day])])$out) != 0){
      orthomosaic <- orthomosaic[-which(orthomosaic$tree_id == tree &
                                          orthomosaic$date == days[day] &
                                          orthomosaic$ndvi %in% c(boxplot.stats(orthomosaic$ndvi[which(orthomosaic$tree_id == tree &
                                                                                                         orthomosaic$date == days[day])])$out)),]
    }
  }
}

# save resulting data frame
save(orthomosaic, file = "out/orthomosaic/outlier_free_ndvi_long_format_phenoclasses_orthomosaic.RData")

# after removing outliers, calculate daily mean NDVI values
orthomosaic_mean_per_tree <- orthomosaic %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

# save resulting data frame
save(orthomosaic_mean_per_tree, file = "out/orthomosaic/outlier_free_ndvi_daily_mean_per_tree_orthomosaic.RData")

## TreeTalker
# remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(treetalker$tree_id)){
  days <- unique(treetalker$date[which(treetalker$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(treetalker$ndvi[which(treetalker$tree_id == tree & 
                                                  treetalker$date == days[day])])$out) != 0){
      treetalker <- treetalker[-which(treetalker$tree_id == tree & 
                                        treetalker$date == days[day] & 
                                        treetalker$ndvi %in% c(boxplot.stats(treetalker$ndvi[which(treetalker$tree_id == tree & 
                                                                                                     treetalker$date == days[day])])$out)),]
    }
  }
}

# after removing outliers, calculate daily mean NDVI values
treetalker_daily_mean <- treetalker %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

# save resulting data frame
save(treetalker_daily_mean, file = "out/treetalker/outlier_free_ndvi_daily_mean_per_tree_treetalker.RData")


#############################################
## create different all-in-one dataframes ##
###########################################
# ----------------------------------------------------------------------------
## all-in-one daily mean NDVI per tree ---------------------------------------

tt_long_phenoclasses %>% relocate(timestamp, .after = date)

load("out/planetscope/outlier_free_daily_ndvi_mean_per_tree_planetscope.RData");planetscope_mean <- as.data.frame(ungroup(planetscope_mean_per_tree));rm(planetscope_mean_per_tree)
load("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.RData");sentinel2_mean <- rename(ndvi_long_pheno_sentinel2, ndvi_mean = ndvi);rm(ndvi_long_pheno_sentinel2);sentinel2_mean$ndvi_sd <- rep(NA, nrow(sentinel2_mean));sentinel2_mean <- sentinel2_mean %>% relocate(ndvi_sd, .after = ndvi_mean)
load("out/treetalker/outlier_free_ndvi_daily_mean_per_tree_treetalker.RData");treetalker_mean <- as.data.frame(ungroup(treetalker_daily_mean));rm(treetalker_daily_mean)
load("out/orthomosaic/outlier_free_ndvi_daily_mean_per_tree_orthomosaic.RData");orthomosaic_mean <- as.data.frame(ungroup(orthomosaic_mean_per_tree));rm(orthomosaic_mean_per_tree)

planetscope_mean$platform <- rep("planetscope", nrow(planetscope_mean))
sentinel2_mean$platform <- rep("sentinel2", nrow(sentinel2_mean))
treetalker_mean$platform <- rep("treetalker", nrow(treetalker_mean))
orthomosaic_mean$platform <- rep("orthomosaic", nrow(orthomosaic_mean))

aio_daily_ndvi_means <- rbind(planetscope_mean, sentinel2_mean, treetalker_mean, orthomosaic_mean)

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

# =============================================================
