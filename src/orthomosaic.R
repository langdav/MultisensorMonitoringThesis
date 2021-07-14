library(raster);library(rgeos);library(rgdal); library(sf); library(RStoolbox)

####################################################################################
## Micasense Channels: Blue, Green, Red, Red Edge, Near Infrared, FLIR (thermal) ##
##################################################################################
# read in 16bit raster (value range 0 - 65535 instead of the value range of 0 - 255 of a 8bit image)
test <- stack("data/orthomosaic/2021_05_14_orthomosaic_new-2-1.tif")
plotRGB(test, r = 3, g = 2, b = 1)
plot(test$X2021_05_14_orthomosaic_new.0.2.1)

#crop to smaller (test) extent
test_shp <- read_sf("data/orthomosaic/test.gpkg")
test_shp <- st_transform(test_shp, 25832)
test_small <- crop(test, test_shp)

# save 16 bit cropped raster
raster::writeRaster(test_small, "data/orthomosaic/test_small.tif", overwrite = T)

# create and save 8bit raster
test_small_8b <- raster::stretch(test_small, minv = 0, maxv = 510,filename = 'data/orthomosaic/test_small_8b.tif', datatype='INT1U')

# read in cropped image
test_small <- stack("data/orthomosaic/test_small.tif")

#convert to 0-255 using the calc. function and basic raster algebra (takes relatively long)
test_small_8b2 <- calc(test_small, fun=function(x){((x - min(x)) * 255)/(max(x)- min(x)) + 0})

######################################################
## merge multiple raster tiles into one big raster ##
####################################################
# files <- list.files("data/orthomosaic/2021_05_04/")
# x <- list()
# for(i in 1:length(files)){
#   tmp <- stack(paste0("data/orthomosaic/2021_05_04/", files[i]))
#   x[[i]] <- tmp
#   rm(tmp)
# }
# m <- do.call(merge, x) #takes ages to process

######################################
## calculate NDVI for single trees ##
####################################
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
trees <- st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

## calculate ndvi for all trees
ndvi_all_trees <- list()
for(i in 1:nrow(trees)){
  for(u in 1:length(x)){
    if(!is.na(raster::cellFromXY(x[[u]], xy = c(trees$easting[i],trees$northing[i])))){
      las_shp <- rgdal::readOGR(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg"))
      crs(las_shp) <- CRS("+proj=longlat +datum=WGS84")
      single_tree <- crop(x[[u]], las_shp)
      ndvi <- RStoolbox::spectralIndices(single_tree, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
      ndvi_all_trees[[i]] <- ndvi
    }
  }
}
#saveRDS(ndvi_all_trees, paste0("data/orthomosaic/ndvi/ndvi_per_tree_20210504.rds"))
#plot(ndvi_all_trees[[24]])






################################################################
## Blue, Green, Red, Red Edge, Near Infrared, FLIR (thermal) ##
##############################################################
plot(test_small$test_small.1)

# calculate NDVI
ndvi <- RStoolbox::spectralIndices(test, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
plot(ndvi)

# crop to single tree extent
load("data/trees_all.RData")
trees <- st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
#trees <- trees %>% filter(tree_id %in% c("mof_cst_00001","mof_cst_00003","mof_cst_00006","mof_cst_00013","mof_cst_00032","mof_cst_00036","mof_cst_00050", "BSF_1"))
trees <- trees[c(1:50),] #reduce to trees with a budburst record

las_shp <- rgdal::readOGR(paste0("data/single_tree_shapefiles/",trees$tree_id[8],".gpkg"))
crs(las_shp) <- CRS("+proj=longlat +datum=WGS84")

plot(test_small$test_small.1); plot(las_shp, add = T)


test_small_crop <- crop(test_small, las_shp)

plotRGB(test_small_crop, r = "test_small.1", g = "test_small.2", b = "test_small.3", scale = 32767)


ndvi <- RStoolbox::spectralIndices(test_small_crop, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
plot(ndvi)
