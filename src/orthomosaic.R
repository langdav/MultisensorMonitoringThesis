#author: David Langenohl
#last modified: 17.08.2021
#description: calculate NDVI values from Orthomosaics
#NOTE: Micasense Channels in Order: Blue, Green, Red, Red Edge, Near Infrared, FLIR (thermal)

rm(list = ls())
library(raster);library(rgeos);library(rgdal); library(sf); library(RStoolbox)

## list available mosaics
mosaics <- list.files("data/orthomosaic/", pattern = ".tif")

## calculate ndvi for all trees
#ndvi_all_trees <- list()
ndvi_long <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)
ndvi_mean_per_tree <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)

for(mosaic in mosaics){
  ortho <- stack(paste0("data/orthomosaic/", mosaic))
  
  load("data/trees.RData")
  trees <- st_transform(trees, 25832)
  trees$tree_id <- as.character(trees$tree_id)
  trees <- trees[c(1:50),] #reduce to trees with a budburst record
  crs(ortho) <- crs(trees)
  
  for(i in 1:nrow(trees)){
    las_shp <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg")) #load single tree shapefile
    single_tree <- crop(ortho, las_shp) #crop tile to single tree extent
    ndvi <- RStoolbox::spectralIndices(single_tree, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
    #ndvi_all_trees[[i]] <- ndvi #write single tree NDVI-rasters into list
    ndvi <- ndvi@data@values
    
    # as in some mosaics not all trees lay within the extent, trees that contain more than 10 percent NA get removed
    if(sum(is.na(ndvi)) < 0.1*length(ndvi)){
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
}


## merge with phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

ndvi_long_pheno_orthomosaic <- merge(ndvi_long, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F)
ndvi_mean_pheno_orthomosaic <- merge(ndvi_mean_per_tree, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F)

## save 
#write.csv(ndvi_long_pheno, "out/orthomosaic/ndvi_long_format_phenoclasses_orthomosaic.csv", row.names = F) #too big for GitHub; save as RData instead
save(ndvi_long_pheno_orthomosaic, file = "out/orthomosaic/ndvi_long_format_phenoclasses_orthomosaic.RData")
save(ndvi_mean_pheno_orthomosaic, file = "out/orthomosaic/ndvi_mean_per_tree_long_format_phenoclasses_orthomosaic.RData")
