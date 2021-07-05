library(dplyr);library(ggplot2);library(viridis);library(tidyr)
library(raster);library(rgeos);library(rgdal)

## read in budburst data
budburst <- read.csv("data/budburst_data/budburst_long.csv")
budburst$date <- as.Date(budburst$date)
budburst$phenoclass <- factor(budburst$phenoclass, levels = c("no_budburst","budburst"))

## NDVI from planetscope data 
# turn off factors
options(stringsAsFactors = FALSE)

# import the planetscope data
# 1 - Blue, 2 - Green, 3 - Red, 4 - NIR
ps_data_st <- stack("data/satellite_data/planetscope/20210307_102818_56_2416_3B_AnalyticMS_SR_clip.tif")

# NDVI = (B4-B3)/(B4+B3)
ndvi <- (ps_data_st[[4]] - ps_data_st[[3]])/(ps_data_st[[4]] + ps_data_st[[3]])
plot(ndvi)


