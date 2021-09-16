#author: David Langenohl; based on the work of Dominic Fawcett
#last modified: 19.08.2021
#description: checking for correlations between PlanetScope and TreeTalker NDVI

rm(list = ls())
#library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);library(minpack.lm) #minpack.lm isneeded for the use of minpack.lm::nlsLM()
library(ggplot2)

#load data
load("out/all_in_one/aio_daily_ndvi_all.RData")
aio_all <- aio_daily_ndvi_all;rm(aio_daily_ndvi_all)
aio_all <- aio_all[-which(aio_all$date > as.Date("2021-07-01")),] #limit to data before "2021-07-01"
aio_all <- aio_all[which(aio_all$platform %in% c("planetscope","treetalker")),] #limit to planetscope and treetalker data

load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
aio_mean <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
aio_mean <- aio_mean[-which(aio_mean$date > as.Date("2021-07-01")),] #limit to data before "2021-07-01"
aio_mean <- aio_mean[which(aio_mean$platform %in% c("planetscope","treetalker")),] #limit to planetscope and treetalker data

load("out/all_in_one/aio_daily_ndvi_per_tree_medians.RData")
aio_median <- aio_daily_ndvi_medians;rm(aio_daily_ndvi_medians)
aio_median <- aio_median[-which(aio_median$date > as.Date("2021-07-01")),] #limit to data before "2021-07-01"
aio_median <- aio_median[which(aio_median$platform %in% c("planetscope","treetalker")),] #limit to planetscope and treetalker data

# plotting; mean NDVI
ggplot(aio_mean[which(aio_mean$tree_id %in% unique(aio_mean$tree_id)[1:5]),], aes(x=date, y=ndvi_mean, colour=platform)) +
  geom_errorbar(aes(ymin=ndvi_mean-ndvi_sd, ymax=ndvi_mean+ndvi_sd), width=.1) +
  geom_point() +
  geom_smooth(formula = y ~ x) +
  facet_grid(tree_id~., scales = "free_y")






#############################
## Spearman correlation coefficient
cor.test(aio_mean$ndvi_mean[which(aio_mean$platform=="planetscope")], aio_mean$ndvi_mean[which(aio_mean$platform=="treetalker")], method = "spearman") # p < 0.05 -> correlation; rho > 0.5 -> große Effektstärke
cor.test(plan_ndvi$budburst_perc, plan_ndvi$ndvi, method = "spearman", alternative = "greater") # p < 0.05 -> positive correlation
cor.test(plan_ndvi$budburst_perc, plan_ndvi$ndvi, method = "spearman", alternative = "less") # p = 1 -> not a negative correlation