rm(list=ls())
library(sf)
library(raster)
library(lidR)
#if(!require("lidRPlugins")) remotes::install_github("Jean-Romain/lidRplugins")
library(lidRplugins) #for additional lidR functions
library(itcSegment)
library(alphashape3d)
library(spatstat)
#if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#if(!require("EBImage")) BiocManager::install("EBImage")
library(EBImage) #needed for watershed
library(RSDB)
#if(!require("lidRPlugins")) remotes::install_github("hunzikp/velox")
#library(velox)

las_all <- readLAS("data/densecloud/2021_03_23_densecloud.laz")
proj4string(las_all) <- CRS("+init=epsg:25832")
#plot(las_all)
plot(las_all, backend = "lidRviewer") 

test_shp <- read_sf("data/orthomosaic/test.gpkg")
test_shp <- st_transform(test_shp, 25832)
las_small <- clip_roi(las_all, test_shp)

  
trees <- readRDS("data/trees.RDS")
tree <- trees[3,]

#create extent to clip las to an area of 70m around the tree of interest
ext <- RSDB::extent_diameter(st_coordinates(tree)[1,1],
                             st_coordinates(tree)[1,2],
                             70)

#clip las to extent
las_small <- clip_roi(las_all, ext)
plot(las_small) # detection of single trees more or less useless here; the point clouds only show the crowns upper surface, nothing else
