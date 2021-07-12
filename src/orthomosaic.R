library(raster);library(rgeos);library(rgdal)

test <- stack("data/orthomosaic/2021_05_14_orthomosaic.tif")
## R G B RED NIR FLIR 
## alternativ: Oder B G R RED NIR FLIR
plot(test$X2021_05_14_orthomosaic.4)
