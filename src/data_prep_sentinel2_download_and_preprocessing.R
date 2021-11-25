#author: David Langenohl
#original script by: Nico Friess (friess@staff.uni-marburg.de)
#last modified: 25.11.2021
#description: Query Copernicus hub, download sentinel scenes and preprocess
rm(list=ls())

## load packages
#remotes::install_github("envima/Rsenal2")
library(Rsenal2);library(tidyverse);
library(reticulate)

## prepare python and import sentinelsat module
if(!reticulate:::python_has_module(python = reticulate::py_config()[[1]],module = "sentinelsat")){
  stop("You need to install the python module 'sentinelsat' on your system.")
}

## initiate Sentinelsat-API with Copernicus credentials
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

## Query for available sentinel 2 scenes ####
product1 = api$query(area = footprint, date = c("20210301","20210531"),
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


