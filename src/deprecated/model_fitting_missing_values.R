#author: David Langenohl
#last modified: 17.08.2021
#description: 
#NOTE: 
rm(list = ls())
library(raster);library(rgeos);library(rgdal);library(sf);library(RStoolbox);library(plyr);library(dplyr)

#load data
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/median_fitted_models_output.RData")

load("data/trees.RData")
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

#find missing SOS values
mean_missing <- model_fitting_out_mean[which(is.na(model_fitting_out_mean$SOS)),]
median_missing <- model_fitting_out_median[which(is.na(model_fitting_out_median$SOS)),]

#merge by tree_id and platform; which missing values are the same in mean and median
both_missing <- merge(mean_missing[,c(1:5,10)], median_missing[,c(1:5)], by = c("tree_id", "platform"), all = F)

unique(both_missing$tree_id)
length(unique(both_missing$tree_id))
setdiff(trees$tree_id,unique(both_missing$tree_id)) #which trees do have SOS values in every platform

#find out, if any trees are not laying within an images extent and are, therefore, NA (which would then be the reason for missing SOS)
#ortho
mosaics <- list.files("data/orthomosaic/", pattern = ".tif")
for(mosaic in mosaics){
  ortho <- stack(paste0("data/orthomosaic/", mosaic))
  for(i in unique(both_missing$tree_id)){
    single_tree_sf <- sf::read_sf(paste0("data/single_tree_shapefiles/",i,".gpkg"))
    single_tree_sf <- single_tree_sf$geom %>% 
      st_union() %>% 
      st_convex_hull()
    single_tree_sf <- st_as_sf(single_tree_sf)
    green <- raster::extract(ortho[[2]], single_tree_sf, na.rm = F)[[1]]
    print(paste0(mosaic,", tree ",i,": number of NAs = ",length(which(is.na(green)))))
  }
}

#Sentinel2
files_sentinel2 <- list.files("data/satellite_data/sentinel2/",pattern = ".tif")
red_stack <- stack(paste0("data/satellite_data/sentinel2/", files_sentinel2[9])) 
red_stack <- raster::dropLayer(red_stack, which(substr(names(red_stack),2,5) == "2020")) #drop layers from 2020
red_stack <- raster::dropLayer(red_stack, 1) #drop layers before the 15.03.

for(u in 1:length(names(red_stack))){
  for(i in unique(both_missing$tree_id)){
    single_tree_sf <- sf::read_sf(paste0("data/single_tree_shapefiles/",i,".gpkg"))
    single_tree_sf <- sf::st_transform(single_tree_sf ,st_crs(red_stack)$proj4string)
    single_tree_sf <- single_tree_sf$geom %>% 
      st_union() %>% 
      st_convex_hull()
    single_tree_sf <- st_as_sf(single_tree_sf)
    green <- raster::extract(red_stack[[u]], single_tree_sf, na.rm = F)[[1]]
    print(paste0(names(red_stack)[u],", tree ",i,": number of NAs = ",length(which(is.na(green)))))
    
  }
}
