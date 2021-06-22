library(raster)

files <- list.files("data/raw/Sentinel1")

test <- raster::raster(paste0("data/raw/Sentinel1/", files[4]))
crs(test) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
plot(test)

mof <- sf::st_read("data/mof_extent_wgs84.gpkg")
#mof <- sf::st_transform(mof, crs(test))
#sf::st_crs(mof)
#sf::write_sf(mof, "data/mof_extent_test.gpkg")

test_small <- crop(test, extent(mof))
plot(test_small)

date <- substring(files[1], 15, 22)

raster::writeRaster(test_small, paste0("data/raw/Sentinel1_cut/sentinel1_", date, ".tiff"), overwrite = T)


cr <- crop(test, mof, snap="out")                    
fr <- rasterize(mof, cr)   
lr <- mask(x=cr, mask=fr)

raster::writeRaster(lr, paste0("data/raw/Sentinel1_cut/sentinel1_", date, ".tiff"), overwrite = T)


leck_mich <- raster(paste0("data/raw/Sentinel1_cut/sentinel1_", date, "_manual.tif"))
plot(leck_mich)
