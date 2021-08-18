#author: David Langenohl
#last modified: 17.08.2021
#description: calculate NDVI values from Orthomosaics
#NOTE: Micasense Channels in Order: Blue, Green, Red, Red Edge, Near Infrared, FLIR (thermal)

rm(list = ls())
library(raster);library(rgeos);library(rgdal);library(sf);library(RStoolbox);library(plyr);library(dplyr)

## list available mosaics
mosaics <- list.files("data/orthomosaic/", pattern = ".tif")

## calculate ndvi for all trees
ndvi_all <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)

for(mosaic in mosaics){
  ortho <- stack(paste0("data/orthomosaic/", mosaic))
  
  load("data/trees.RData")
  trees <- st_transform(trees, 25832)
  trees$tree_id <- as.character(trees$tree_id)
  trees <- trees[c(1:50),] #reduce to trees with a budburst record
  crs(ortho) <- crs(trees)
  
  for(i in 1:nrow(trees)){
    cat("Processing", trees_tree_id[i], "in Orthomosaic", mosaic,"\n")
    
    #load single tree shapefile
    single_tree_sf <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg"))
    
    #extract green for single tree; do not remove NAs to check if tree lies within the boundaries of the orthomosaic
    green <- raster::extract(ortho[[2]], single_tree_sf, na.rm = F)
    
  
    # as in some mosaics not all trees lay within the extent, trees that contain more than 10 percent NA get removed
    if(sum(is.na(green)) < 0.1*length(green)){
      #remove NAs from green (if present)
      green <- na.omit(green)
      
      #extract  red and NIR values for single tree 
      red <- raster::extract(ortho[[3]], single_tree_sf, na.rm = T)
      nir <- raster::extract(ortho[[5]], single_tree_sf, na.rm = T)
      
      #calculate NDVI = (nir-red)/(nir+red)
      ndvi <- (nir-red)/(nir+red)
      
      #compute 80th percentile of "green" values
      green80thpercentile <- as.numeric(quantile(green, 0.8, na.rm=T))
      
      #perform brightest pixel masking of NDVI values for reduced background effect; remove NDVI values of cells, where the green band is smaller than the 80th percentile
      ndvi_80 <- ndvi[which(green>green80thpercentile)]
      
      #create output data frame
      ndvi_all <- rbind(ndvi_all, data.frame(tree_id=rep(trees$tree_id[i], length(ndvi_80)),
                                               date=rep(as.Date(substr(mosaic, 1, 10), "%Y_%m_%d"), length(ndvi_80)),
                                               ndvi = ndvi_80))
    }
  }
}

#merge with phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

ndvi_all_pheno <- merge(ndvi_all, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F)

#save resulting data frame
save(ndvi_all_pheno, file = "out/orthomosaic/ndvi_all_with_phenoclasses_orthomosaic.RData")

#remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(ndvi_all_pheno$tree_id)){
  cat("Processing", tree, "\n")
  days <- unique(ndvi_all_pheno$date[which(ndvi_all_pheno$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(ndvi_all_pheno$ndvi[which(ndvi_all_pheno$tree_id == tree &
                                                   ndvi_all_pheno$date == days[day])])$out) != 0){
      ndvi_all_pheno <- ndvi_all_pheno[-which(ndvi_all_pheno$tree_id == tree &
                                          ndvi_all_pheno$date == days[day] &
                                          ndvi_all_pheno$ndvi %in% c(boxplot.stats(ndvi_all_pheno$ndvi[which(ndvi_all_pheno$tree_id == tree &
                                                                                                         ndvi_all_pheno$date == days[day])])$out)),]
    }
  }
}

# save resulting data frame
save(ndvi_all_pheno, file = "out/orthomosaic/outlier_free_ndvi_all_with_phenoclasses_orthomosaic.RData")

# after removing outliers, calculate daily mean NDVI values
ndvi_per_tree <- ndvi_all_pheno %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

# save resulting data frame
save(ndvi_per_tree, file = "out/orthomosaic/outlier_free_ndvi_mean_per_tree_with_phenoclasses_orthomosaic.RData")

