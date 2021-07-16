library(plyr);library(dplyr);library(ggplot2);library(viridis);library(tidyr)
library(raster);library(rgeos);library(rgdal);library(lidR);library(sf)

## load data
# trees <- readRDS("data/trees.RDS")
load("data/trees_all.RData")
trees <- st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
#trees <- trees %>% filter(tree_id %in% c("mof_cst_00001","mof_cst_00003","mof_cst_00006","mof_cst_00013","mof_cst_00032","mof_cst_00036","mof_cst_00050", "BSF_1"))
trees <- trees[c(1:50),] #reduce to trees with a budburst record

## create colorblind friendle palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

############################################################
## create NDVI (or other indices) from Sentinel-2 images ##
##########################################################
## load mof extent shp file to crop sentinel image
las_shp <- sf::read_sf("data/mof_extent/mof_extent_wgs84.gpkg")
las_shp <- sf::st_transform(las_shp,  "+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") #reproject to EPSG 42106

files_sentinel2 <- list.files("data/satellite_data/sentinel2/",pattern = ".tif")
substr(files_sentinel2, 32, 34)

## load and crop red channel
red <- stack(paste0("data/satellite_data/sentinel2/", files_sentinel2[9])) #import data; EPSG: 42106
red <- raster::dropLayer(red, which(substr(names(red),2,5) == "2020")) #drop layers from 2020
red <- raster::dropLayer(red, c(1,9:13)) #drop layers beyond the 15.03.-01.06. time period
red <- raster::dropLayer(red, c(3:5,7)) #drop layers without data, due to clouds
red <- crop(red, las_shp)

## load and crop nir channel
nir <- stack(paste0("data/satellite_data/sentinel2/", files_sentinel2[5])) #import data; EPSG: 42106
nir <- raster::dropLayer(nir, which(substr(names(nir),2,5) == "2020")) #drop layers from 2020
nir <- raster::dropLayer(nir, c(1,9:13))  #drop layers beyond the 15.03.-01.06. time period
nir <- raster::dropLayer(nir, c(3:5,7)) #drop layers without data, due to clouds
nir <- crop(nir, las_shp)

## NDVI of whole mof
#NDVI = (nir-red)/(nir+red)
ndvi <- (nir-red)/(nir+red)
plot(ndvi)

## single tree NDVI; raster
ndvi_all_trees <- list()
for(i in 1:nrow(trees)){
  las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg")) #load single tree shapefile
  sf::st_crs(las_shp) <- 25832
  las_shp <- sf::st_transform(las_shp,  "+proj=laea +lat_0=55 +lon_0=20 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") #reproject to EPSG 42106
  single_tree_red <- crop(red, las_shp) #crop tile to single tree extent
  single_tree_nir <- crop(nir, las_shp) #crop tile to single tree extent
  ndvi <- (single_tree_nir-single_tree_red)/(single_tree_nir+single_tree_red)
  ndvi_all_trees[[i]] <- ndvi #write single tree NDVIs into list
  
}
#saveRDS(ndvi_all_trees, "out/sentinel2/ndvi_per_tree_raster.rds")
#plot(ndvi_all_trees[[24]])

## single tree NDVI; values
dates <- c(as.Date("2021-03-20"),
           as.Date("2021-03-30"),
           as.Date("2021-05-09"))
ndvi_all_trees_long <- data.frame()

for(i in 1:length(ndvi_all_trees)){
  values <- c(mean(ndvi_all_trees[[i]]@data@values[,1]),
              mean(ndvi_all_trees[[i]]@data@values[,2]),
              mean(ndvi_all_trees[[i]]@data@values[,3]))
  tmp_ndvi_df <- data.frame(tree_id = rep(trees$tree_id[i], 3),
                      date = dates,
                      ndvi = values)
  ndvi_all_trees_long <- rbind(ndvi_all_trees_long, tmp_ndvi_df)
}

## merge with phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

ndvi_sentinel2_phenoclasses <- merge(ndvi_all_trees_long, phenoclasses, by = c("tree_id","date"), all.x = F, all.y = F) #merge with phenoclasses
#write.csv(ndvi_sentinel2_phenoclasses, "out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.csv")

## plotting
ndvi_sentinel2_phenoclasses %>% 
  ggplot(aes(x=date, y=ndvi, group=tree_id, color=as.factor(budburst_perc))) +
  geom_point(size = 3)+
  geom_line() +
  #facet_grid(tree_id~., scales = "free_y") +
  #scale_color_manual(values= cbPalette) +
  scale_color_viridis(discrete = T) +
  scale_x_date(date_breaks = "1 week") +
  xlab("Date") +
  ylab("NDVI") +
  labs(color="buds bursted (in %)") +
  theme_light()

