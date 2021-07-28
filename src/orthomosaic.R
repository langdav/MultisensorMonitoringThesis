rm(list = ls())
library(raster);library(rgeos);library(rgdal); library(sf); library(RStoolbox)

####################################################################################
## Micasense Channels: Blue, Green, Red, Red Edge, Near Infrared, FLIR (thermal) ##
##################################################################################
# read in 16bit raster (value range 0 - 65535 instead of the value range of 0 - 255 of a 8bit image)
test <- stack("data/orthomosaic/2021_04_23_orthomosaic.tif")

load("data/trees.RData")
trees <- sf::st_transform(trees, 25832)
crs(test) <- crs(trees)

plot(test$X2021_04_23_orthomosaic.1);plot(trees$geometry, add = T, col = "red", lwd = 15)


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
load("data/trees.RData")
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
load("data/trees.RData")
trees <- st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record
trees <- trees[-c(30,37),] #those two trees lay beyond the mosaics boundaries and, hence, need to be excluded

## list available mosaics
mosaics <- list.files("data/orthomosaic/", pattern = ".tif")
  
## calculate ndvi for all trees
#ndvi_all_trees <- list()
ndvi_long <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)
ndvi_mean_per_tree <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)

for(mosaic in mosaics){
  ortho <- stack(paste0("data/orthomosaic/", mosaic))
  for(i in 1:nrow(trees)){
    las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg")) #load single tree shapefile
    # if(is.na(crs(las_shp))){
    #   sf::st_crs(las_shp) <- 25832
    # }
    single_tree <- crop(ortho, las_shp) #crop tile to single tree extent
    ndvi <- RStoolbox::spectralIndices(single_tree, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
    #ndvi_all_trees[[i]] <- ndvi #write single tree NDVI-rasters into list
    ndvi <- ndvi@data@values
    ndvi_long <- rbind(ndvi_long, data.frame(tree_id=rep(trees$tree_id[i], length(ndvi)),
                                             date=rep(as.Date(substr(mosaic, 1, 10), "%Y_%m_%d"), length(ndvi)),
                                             ndvi = ndvi))
    
    ## remove outliers and calculate NDVI mean values per tree
    outliers <- which(ndvi %in% c(boxplot.stats(ndvi)$out)) 
    outliers_removed <- ndvi[-outliers] #remove rows that are containing the outliers
    
    if(length(outliers_removed) != 0){
      ndvi_mean <- mean(outliers_removed)
      } else {
      ndvi_mean <- mean(ndvi)
      }
    
    
    ndvi_mean_per_tree <- rbind(ndvi_mean_per_tree, data.frame(tree_id=trees$tree_id[i],
                                                               date=as.Date(substr(mosaic, 1, 10), "%Y_%m_%d"),
                                                               ndvi = ndvi_mean))
    
  }
}
#saveRDS(ndvi_all_trees, paste0("out/orthomosaic/ndvi_per_tree_20210504.rds"))


## merge with phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

ndvi_long_pheno <- merge(ndvi_long, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F)
ndvi_mean_pheno <- merge(ndvi_mean_per_tree, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F)

## save 
write.csv(ndvi_long_pheno, "out/orthomosaic/ndvi_long_format_phenoclasses_orthomosaic.csv", row.names = F)
write.csv(ndvi_mean_pheno, "out/orthomosaic/ndvi_mean_per_tree_long_format_phenoclasses_orthomosaic.csv", row.names = F)

