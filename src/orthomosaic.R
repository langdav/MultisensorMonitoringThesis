rm(list = ls())
library(raster);library(rgeos);library(rgdal); library(sf); library(RStoolbox)

####################################################################################
## Micasense Channels: Blue, Green, Red, Red Edge, Near Infrared, FLIR (thermal) ##
##################################################################################
# read in 16bit raster (value range 0 - 65535 instead of the value range of 0 - 255 of a 8bit image)
test <- stack("data/orthomosaic/2021_05_14_orthomosaic_new-2-1.tif")
test <- stack("data/orthomosaic/2021_04_23_orthomosaic.tif")
plotRGB(test, r = 3, g = 2, b = 1)
plot(test$X2021_04_23_orthomosaic.4);plot(las_shp, add=T);plot(single_tree, add=T);plot(single_tree2, add=T)

plot(las_shp);plot(single_tree, add=T)


##############################################################
i <- 2
las_shp <- sf::read_sf(paste0(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg"))) #load single tree shapefile)
plot(las_shp$geom);plot(trees$geometry[i], add = T, col = "red", lwd = 15)

plot(trees$geometry[1],  col = "red", lwd = 15)


load("data/trees.RData")
trees <- sf::st_transform(trees, 25832)
crs(test) <- crs(trees)

plot(test);plot(trees$geometry, add = T, col = "red", lwd = 15)

##############################################################

las_shp <- sf::read_sf(paste0("data/mof_extent/mof_extent_new2.gpkg")) #load single tree shapefile
las_shp <- sf::st_transform(las_shp,crs(test))
test2 <- crop(test$X2021_04_27_orthomosaic.4, las_shp) #
plot(test2);plot(las_shp, add=T)

#crop to smaller (test) extent
test_shp <- read_sf("data/orthomosaic/test.gpkg")
test_shp <- st_transform(test_shp, 25832)
test_small <- crop(test, test_shp)
plot(test_small$test.1)

# save 16 bit cropped raster
raster::writeRaster(test_small, "data/orthomosaic/test_small.tif", overwrite = T)

# create and save 8bit raster
test_small_8b <- raster::stretch(test_small, minv = 0, maxv = 510,filename = 'data/orthomosaic/test_small_8b.tif', datatype='INT1U')

# read in cropped image
test_small <- stack("data/orthomosaic/test_small.tif")

#convert to 0-255 using the calc. function and basic raster algebra (takes relatively long)
test_small_8b2 <- calc(test_small, fun=function(x){((x - min(x)) * 255)/(max(x)- min(x)) + 0})

#######################################################################
## merge multiple raster tiles into one big raster; NOT RECOMMENDED ##
#####################################################################
# files <- list.files("data/orthomosaic/2021_05_04/")
# x <- list()
# for(i in 1:length(files)){
#   tmp <- stack(paste0("data/orthomosaic/2021_05_04/", files[i]))
#   x[[i]] <- tmp
#   rm(tmp)
# }
# m <- do.call(merge, x) #takes ages to process

#########################################################
## calculate NDVI for single trees using raster tiles ##
#######################################################
## get tiles into a list
files <- list.files("data/orthomosaic/2021_05_04/")
x <- list()
for(i in 1:length(files)){
  tmp <- stack(paste0("data/orthomosaic/2021_05_04/", files[i]))
  x[[i]] <- tmp
  rm(tmp)
}

## load tree data 
load("data/trees_all.RData")
trees <- sf::st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

## calculate ndvi for all trees
ndvi_all_trees <- list()
for(i in 1:nrow(trees)){
  for(u in 1:length(x)){
    if(!is.na(raster::cellFromXY(x[[u]], xy = c(trees$easting[i],trees$northing[i])))){
      las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg")) #load single tree shapefile
      sf::st_crs(las_shp) <- 25832 #set crs
      las_shp <- sf::st_transform(las_shp,  3857) #transform crs to crs of orthomosaics
      single_tree <- crop(x[[u]], las_shp) #crop tile to single tree extent
      ndvi <- RStoolbox::spectralIndices(single_tree, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
      ndvi_all_trees[[i]] <- ndvi #write single tree NDVIs into list
    }
  }
}
#saveRDS(ndvi_all_trees, paste0("out/orthomosaic/ndvi_per_tree_20210504.rds"))
#plot(ndvi_all_trees[[24]])


#############################################################
## calculate NDVI for single trees using a single BigTIFF ##
###########################################################
## load tree data 
load("data/trees_all.RData")
trees <- st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

## calculate ndvi for all trees
ndvi_all_trees <- list()
for(i in 1:nrow(trees)){
  las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[19],".gpkg")) #load single tree shapefile
  if(is.na(crs(las_shp))){
    sf::st_crs(las_shp) <- 25832
  }
  single_tree <- crop(test, las_shp) #crop tile to single tree extent
  ndvi <- RStoolbox::spectralIndices(single_tree, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
  ndvi_all_trees[[i]] <- ndvi #write single tree NDVIs into list
}
#saveRDS(ndvi_all_trees, paste0("out/orthomosaic/ndvi_per_tree_20210504.rds"))
#plot(ndvi_all_trees[[24]])

#########################################################################
## calculate mean NDVI for single trees and match witch budburst data ##
#######################################################################
# phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

ndvi_mean <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)
for(i in 1:nrow(trees)){
  las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg")) #load single tree shapefile
  if(is.na(crs(las_shp))){
    sf::st_crs(las_shp) <- 25832
  }
  single_tree <- crop(test, las_shp) #crop tile to single tree extent
  ndvi <- RStoolbox::spectralIndices(single_tree, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
  
  ndvi_mean <- rbind(ndvi_mean, data.frame(tree_id=trees$tree_id[i],
                                           date=as.Date("2021-04-23"),
                                           ndvi = mean(ndvi@data@values)))

}

ndvi_mean_long <- merge(ndvi_mean, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F) #merge with phenoclasses
write.csv(ndvi_mean_long, "out/orthomosaic/ndvi_long_format_phenoclasses_orthomosaic", row.names = F)


#############################################
las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[17],".gpkg")) #load single tree shapefile
sf::st_crs(las_shp) <- 25832
#las_shp <- sf::st_transform(las_shp,  3857)
#las_shp <- sf::st_transform(las_shp,  "+proj=merc +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0") #reproject to EPSG 42106

single_tree <- crop(test, las_shp) #crop tile to single tree extent
#ndvi <- RStoolbox::spectralIndices(single_tree, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI

plotRGB(single_tree, r = 3, g = 2, b = 1)
plot(ndvi)
