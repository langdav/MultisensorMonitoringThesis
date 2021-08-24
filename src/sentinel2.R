rm(list = ls())
library(plyr);library(dplyr);library(ggplot2);library(viridis);library(tidyr)
library(raster);library(rgeos);library(rgdal);library(lidR);library(sf)

# list available files
files_sentinel2 <- list.files("data/satellite_data/sentinel2/",pattern = ".tif")

# load red channel
red_stack <- stack(paste0("data/satellite_data/sentinel2/", files_sentinel2[9])) #import data; EPSG: 42106
red_stack <- raster::dropLayer(red_stack, which(substr(names(red_stack),2,5) == "2020")) #drop layers from 2020
red_stack <- raster::dropLayer(red_stack, 1) #drop layers before the 15.03.

# load nir channel
nir_stack <- stack(paste0("data/satellite_data/sentinel2/", files_sentinel2[5])) #import data; EPSG: 42106
nir_stack <- raster::dropLayer(nir_stack, which(substr(names(nir_stack),2,5) == "2020")) #drop layers from 2020
nir_stack <- raster::dropLayer(nir_stack, 1) #drop layers before the 15.03.

# load 
load("data/trees.RData")
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

ndvi_all <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)
for(i in 1:nrow(trees)){
  cat("Processing", trees$tree_id[i], "\n")
  
  #load single tree shapefile
  las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg"))
  las_shp <- sf::st_transform(las_shp ,st_crs(red_stack)$proj4string) #reproject to projection of Sentinel-2 data (+proj=laea +lat_0=55 +lon_0=20)
  
  # convert Spatial Points to Polygon Shapefile
  las_shp <- las_shp$geom %>% 
    st_union() %>% 
    st_convex_hull()
  las_shp <- st_as_sf(las_shp)
  
  #extract red and NIR values for single tree
  red <- raster::extract(red_stack, las_shp, df = T)
  nir <- raster::extract(nir_stack, las_shp, df = T)
  
  #calculate NDVI = (nir-red)/(nir+red)
  ndvi <- (nir-red)/(nir+red)
  
  #create output data frame
  tree_id_tmp <- NULL
  date_tmp <- NULL
  ndvi_tmp <- NULL
  for(u in 2:ncol(ndvi)){
    tree_id_tmp <- c(tree_id_tmp, rep(trees$tree_id[i], nrow(ndvi)))
    date_tmp <-c(date_tmp, rep(format(as.Date(substr(colnames(ndvi)[u],2,9), "%Y%m%d")), nrow(ndvi)))
    ndvi_tmp <- c(ndvi_tmp, ndvi[,u])
  }
  
  ndvi_all <- rbind(ndvi_all, data.frame(tree_id=tree_id_tmp,
                                         date=as.Date(date_tmp),
                                         ndvi = ndvi_tmp))
}

#remove NA-entries; those data is missing due to clouds
ndvi_all_cloudless <- na.omit(ndvi_all)

#merge with phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

sen2_ndvi_all <- merge(ndvi_all_cloudless, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F)

#save resulting data frame
save(sen2_ndvi_all, file = "out/sentinel2/ndvi_all_with_sentinel2_sentinel2.RData")

#calculate daily mean and median (of applicable) NDVI values
sen2_ndvi_mean_per_tree <- sen2_ndvi_all %>% 
  group_by(tree_id, date) %>% 
  dplyr::summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

sen2_ndvi_median_per_tree <- sen2_ndvi_all %>% 
  group_by(tree_id, date) %>% 
  dplyr::summarize(ndvi_median = median(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

# save resulting data frame
save(sen2_ndvi_mean_per_tree, file = "out/sentinel2/ndvi_mean_per_tree_with_phenoclasses_sentinel2.RData")
save(sen2_ndvi_median_per_tree, file = "out/sentinel2/ndvi_median_per_tree_with_phenoclasses_sentinel2.RData")

