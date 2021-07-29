rm(list = ls())
library(plyr);library(dplyr);library(ggplot2)

# compare daily per tree NDVI

################################################################################
## Planetscope
load("out/planetscope/ndvi_long_format_phenoclasses_planetscope.RData")
planetscope <- ndvi_long_pheno_planetscope; rm(ndvi_long_pheno_planetscope)

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
planetscope_daily_mean <- planetscope %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))


################################################################################
## Sentinel-2
load("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.RData")
sentinel2 <- ndvi_long_pheno_sentinel2; rm(ndvi_long_pheno_sentinel2)
sentinel2$date <- as.Date(sentinel2$date)


################################################################################
## TreeTalker
load("out/treetalker/ndvi_long_format_phenoclasses_treetalker.RData")
treetalker <- ndvi_long_format_phenoclasses_treetalker; rm(ndvi_long_format_phenoclasses_treetalker)

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
