# # # # # # # # # # # # # # # # # # # # # #
# # # # # Script by Nicolas Friess# # # # #
# # # # friess@staff.uni-marburg.de # # # #
# # # # # # 2020-09-09 10:29:26 # # # # # #
# # # # # # # # # # # # # # # # # # # # # #

# Query Copernicus hub, download sentinel scenes and preprocess ####

rm(list=ls())

# load packages ####
#remotes::install_github("envima/Rsenal2")

library(Rsenal2);library(tidyverse);
library(reticulate)
# prepare python and import sentinelsat module ####

if(!reticulate:::python_has_module(python = reticulate::py_config()[[1]],module = "sentinelsat")){
  stop("You need to install the python module 'sentinelsat' on your system.")
}

# Initiate Sentinelsat-API with Copernicus credentials ####

## the user data to access the data from RSDB is stored in a project-specific .Renvrion file
## you can provide your own data with: 
### usethis::edit_r_environ(scope = "project")'; this creates or opens a project-specific .Renviron file
### then insert two lines 'COPERNICUS_USER = username' and 'COPERNICUS_PW = your_password'

sentinelsat = reticulate::import_from_path("sentinelsat", path=Sys.getenv("PYTHONPKG_PATH"))
api = sentinelsat$SentinelAPI(Sys.getenv("COPERNICUS_USER"), Sys.getenv("COPERNICUS_PW"), "https://scihub.copernicus.eu/dhus")

# prepare bounding box of AOI ####
mof <- Rsenal2::mofAreas$mof
mof <- sf::st_transform(mof, crs=4326)
footprint <- sf:::prnt.POLYGON(mof$geom[[1]])

# Query copernicus hub ####

# Sentinel 1 ####
## Query for available sentinel 1 scenes ####
product = api$query(area = footprint, date = c("20210512", "20210515"),
                    platformname = "Sentinel-1", producttype="GRD")
save(product,file="data/satellite_data/sentinel1/products_sentinel1.RData")

## Prepare paths ####

gpt_path <-  normalizePath("C:/Program Files/snap/bin/",mustWork = F)
output_path <-  paste0("\'", normalizePath("data/satellite_data/sentinel1/sentinel1_out//",mustWork = T),"\'")

## Download Sentinel 1 scenes ####

# Here I temporarily download the full sentinel-1 scene ( ~1GB) and at first create a Subset of the scene for our study region
# The steps in the graph are:
#   Open file > Subset scene to study region > save the subset 
# Afterwards the scene is removed from harddrive and subsequent calculations can be done with the subset

for(i in 1:length(product)){
  cat("Starting Scene ",i,"of", length(product),"\n")
  k <- product[[i]]
  if(!gsub(".SAFE","_sub.zip.dim",k$filename) %in% dir("data/satellite_data/sentinel1/sentinel1_inp/")){
    if(product[[i]]$orbitdirection!="ASCENDING") next
    api$download(k$uuid, "data/satellite_data/sentinel1/sentinel1_tmp/")
    input_path <- "data/satellite_data/sentinel1/sentinel1_tmp/"
    input_file <- paste0("\'", normalizePath(input_path), "\\",dir(input_path),"\'")
    #input_file <- paste0(getwd(), "/", input_path, dir(input_path))
    output_file <- gsub("sentinel1_tmp","sentinel1_inp",gsub(".zip","_sub.zip",input_file))
    shell(cmd = paste0("cd ",gpt_path," && ","gpt ",
                       paste0("\'", normalizePath("data/satellite_data/sentinel1/sentinel1_graph/Sentinel1_Subset.xml",mustWork = T),"\'"), 
                       " -Pinput=", 
                       input_file,
                       " -Poutput=",
                       output_file))
    
    unlink(dir(input_path,full.names = T), recursive = T)
  }
  cat("Finished scene", i,"of", length(product), "\n\n\n")
}



################################################################################################################
## processing those files that needed to be downloaded manually (most of them, as most were already offline) ##
##############################################################################################################

whats_there <- list.files("data/raw/Sentinel1_downloads/")

for(i in 1:length(whats_there)){
  cat("Starting Scene ",i,"of", length(whats_there),"\n")
  input_path <- "data/raw/Sentinel1_downloads/"
  input_file <- paste0("\'", normalizePath(input_path), "\\", whats_there[i],"\'")
  #input_file <- paste0(getwd(), "/", input_path, dir(input_path))
  output_file <- paste0(normalizePath("data/satellite_data/sentinel1/sentinel1_inp/"), "\\",gsub(".zip","_sub.zip",whats_there[i]))
  shell(cmd = paste0("cd ",gpt_path," && ","gpt ",
                     paste0("\'", normalizePath("data/satellite_data/sentinel1/sentinel1_graph/Sentinel1_Subset.xml",mustWork = T),"\'"), 
                     " -Pinput=", 
                     input_file,
                     " -Poutput=",
                     output_file))
}

################################################################################################################
################################################################################################################
################################################################################################################

## Run preprocessing as GPT graph in SNAP ####

# The steps in the graph are:
#   Open file > Apply Orbit File > Thermal Noise Filter
#   Remove GRD Border Noise > Calibrate and generate Sigma VV and VH Bands > Apply Speckle Filter (Lee Sigma)
#   Terrain Correction > Linear to Log transform > Write out Tiff
# Afterwards the scene is removed from harddrive and subsequent calculations can be done with the subset
input_path <- "data/satellite_data/sentinel1/sentinel1_inp/"
sen1_scenes <- dir("data/satellite_data/sentinel1/sentinel1_inp/")[grepl(".dim",dir("data/satellite_data/sentinel1/sentinel1_inp/"))]
for(i in 1:length(sen1_scenes)){
  #cat("Starting Scene ",i,"of", length(sen1_scenes),"\n")
  input_file <- paste0("\'", normalizePath(input_path), "\\",sen1_scenes[i],"\'")
  if(!file.exists(gsub("_sub.zip.dim","_sigma_vv.tif",paste0("data/satellite_data/sentinel1/sentinel1_out/",sen1_scenes[i])))){
    output_file <- gsub("sentinel1_inp","sentinel1_out\\",paste0(gsub("_sub.zip.dim\'","",input_file),".tif\'"))
    shell(cmd = paste0("cd ",gpt_path," && ","gpt ",
                       paste0("\'", normalizePath("data/satellite_data/sentinel1/sentinel1_graph/Sentinel1_allout_tif.xml",mustWork = T),"\'"), 
                       " -Pinput=", 
                       input_file,
                       " -Poutput_sigma_vv=",
                       gsub(".tif","_sigma_vv.tif",output_file),
                       " -Poutput_gamma_vv=",
                       gsub(".tif","_gamma_vv.tif",output_file),
                       " -Poutput_sigma_vh=",
                       gsub(".tif","_sigma_vh.tif",output_file),
                       " -Poutput_gamma_vh=",
                       gsub(".tif","_gamma_vh.tif",output_file)
    )
    )
    
  }

  cat("Finished scene", i,"of", length(sen1_scenes), "\n\n\n")
}


# Sentinel 2 ####
# The procedure for sentinel 2 is basically the same as for sentinel 1

## Query for available sentinel 2 scenes ####
product1 = api$query(area = footprint, date = c("20200101","20201023"),
                    producttype = "S2MSI2A", cloudcoverpercentage = c(0,60))
product = api$query(area = footprint, date = c("20200101","20201023"),
                    producttype = "S2MSI2A", cloudcoverpercentage = c(0,90))
product <- product[!names(product)%in%names(product1)]
save(product,file="data/satellite_data/products_sentinel2.RData")

## Prepare paths ####
gpt_path <-  normalizePath("C:/\"Program Files\"/snap/bin",mustWork = F)
graph <- paste0("\'", normalizePath("data/satellite_data/sen_graph/Sentinel2_graph_subout_tif.xml",mustWork = T),"\'")

## Download Sentinel 2 scenes and run GPT graph in SNAP ####

# The GPT graph for sentinel 2 only opens, subsets and writes out the subset, since we already are working with level 2 data (i.e. calibrated and so on)
input_path <- "data/satellite_data/sen_tmp/"
for(i in 1:length(product)){
  k <- product[[i]]
  api$download(k$uuid, "data/satellite_data/sen_tmp")
  input_file <- paste0("\'", normalizePath(input_path),dir(input_path),"\'")
  output_file <- gsub("sen_tmp","sen2_out",paste0(gsub(".zip\'","",input_file),"_sub\'"))
  
  shell(cmd = paste0("cd ",gpt_path," && ","gpt ", graph, " -Pinput=", input_file, " -Poutput=",output_file))
  
  file.remove(dir("data/satellite_data/sen_tmp/",full.names = T))
  cat("Finished scene", i,"of", length(product), "\n\n\n")
}


