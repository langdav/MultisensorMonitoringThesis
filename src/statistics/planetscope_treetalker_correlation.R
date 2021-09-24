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
  geom_errorbar(aes(ymin=ndvi_mean-ndvi_sd, ymax=ndvi_mean+ndvi_sd), width=1) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(tree_id~., scales = "free_y")

# testing for correlations
# Spearman correlation coefficient; insert NAs
plan <- aio_mean[which(aio_mean$platform=="planetscope"),c("tree_id", "date", "ndvi_mean","platform")]
tt <- aio_mean[which(aio_mean$platform=="treetalker"),c("tree_id", "date", "ndvi_mean","platform")]
test <- merge(plan, tt, all = T, by = c("tree_id", "date")) #insert NA for days that have data for one platform but not for the other one

cor.test(test$ndvi_mean.x, test$ndvi_mean.y, method = "spearman") # p < 0.05 -> correlation; rho > 0.5 -> große Effektstärke
cor.test(test$ndvi_mean.x, test$ndvi_mean.y, method = "spearman", alternative = "greater") # p < 0.05 -> positive correlation
cor.test(test$ndvi_mean.x, test$ndvi_mean.y, method = "spearman", alternative = "less") # p = 1 -> not a negative correlation

# Spearman correlation coefficient; predict missing data
library(lubridate)
plan$date_doy <- lubridate::yday(plan$date)
plan_mean <- NULL
for(i in unique(plan$date_doy)){
  ndvi <- mean(plan$ndvi_mean[which(plan$date_doy == i)], na.rm = T)
  plan_mean <- rbind(plan_mean, data.frame(ndvi_mean = ndvi,
                                 date_doy = i))
}

tt$date_doy <- lubridate::yday(tt$date)
tt_mean <- NULL
for(i in unique(tt$date_doy)){
  ndvi <- mean(plan$ndvi_mean[which(plan$date_doy == i)], na.rm = T)
  tt_mean <- rbind(tt_mean, data.frame(ndvi_mean = ndvi,
                                           date_doy = i))
}

plan_model <- lm(ndvi_mean ~ date_doy, plan_mean)
tt_model <- lm(ndvi_mean ~ date_doy, tt_mean)

dates <- seq(lubridate::yday("2021-03-07"),lubridate::yday("2021-06-01"),1)
plan_pred <- as.numeric(predict(plan_model,data.frame(date_doy=dates)))
tt_pred <- as.numeric(predict(tt_model,data.frame(date_doy=dates)))

cor.test(plan_pred, tt_pred, method = "spearman") # p < 0.05 -> correlation; rho > 0.5 -> große Effektstärke
cor.test(plan_pred, tt_pred, method = "spearman", alternative = "greater") # p < 0.05 -> positive correlation
cor.test(plan_pred, tt_pred, method = "spearman", alternative = "less") # p = 1 -> not a negative correlation
