library(sf);library(raster)

sen1_results <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices.RDS")
sen1_pred <- readRDS("data/satellite_data/satellite_results/prediction/sentinel1_predictors.RDS")

plot(sen1_pred[[1]][[2]])
