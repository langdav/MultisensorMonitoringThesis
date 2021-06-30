# # # # # # # # # # # # # # # # # # # # # #
# # # # # Script by Nicolas Friess# # # # #
# # # # friess@staff.uni-marburg.de # # # #
# # # # # # 2020-09-14 10:01:59 # # # # # #
# # # # # # # # # # # # # # # # # # # # # #

rm(list=ls())

# load packages ####
library(tidyverse);library(RSQLite)
library(sf);library(raster);library(RStoolbox);library(mapview)
library(glcm);library(parallel)
# input data ####

#load("data/spatial/trees.RData")
trees <- readRDS("data/trees.RDS")
colnames(trees)[1] <- "tree_id"
load("data/satellite_data/sentinel1/products_sentinel1.RData")
#product[[44]] <- NULL #remove "S1A_IW_GRDH_1SDV_20210331T053413_20210331T053438_037238_0462D6_76C3", as no orbit file available
#product[[23]] <- NULL #remove ""S1B_IW_GRDH_1SDV_20210430T053317_20210430T053342_026692_033029_0BB4", as no orbit file available

# Index calculations ####

dirs <- dir("data/satellite_data/sentinel1/sentinel1_out/",pattern = ".tif" ,full.names = T)
dirs <- split(dirs, as.factor(unlist(lapply(strsplit(basename(dirs),"_"), function(x) x[5]))))
results <- as.list(rep(NA,length(dirs)))
metadata <- unlist(lapply(product, function(x) unlist(strsplit(x$filename,"_"))[5]))
metadata[!metadata%in%names(dirs)]
dirs <- dirs[names(dirs)%in% metadata]

results <- mclapply(as.list(1:length(dirs)), function(i){
  path <- dirs[[i]][1]
  meta <- product[[which(unlist(lapply(product, function(x){isTRUE(all.equal(unlist(strsplit(basename(path),"_"))[1:9] ,unlist(strsplit(gsub(".SAFE","",x$filename),"_"))))})))]]
  dat <- data.frame(tree_id=trees$tree_id,date=lubridate::as_datetime(names(dirs)[i]), orbit = meta$orbitdirection, scene = meta$uuid)
  dat$gamma_vv <- raster::extract(raster(grep("gamma_vv", dirs[[i]],value = T)), sf::st_coordinates(trees))
  dat$gamma_vh <- raster::extract(raster(grep("gamma_vh", dirs[[i]],value = T)), sf::st_coordinates(trees))
  dat$sigma_vv <- raster::extract(raster(grep("sigma_vv", dirs[[i]],value = T)), sf::st_coordinates(trees))
  dat$sigma_vh <- raster::extract(raster(grep("sigma_vh", dirs[[i]],value = T)), sf::st_coordinates(trees))
  
  ind_stack <- raster::stack(raster(grep("gamma_vv", dirs[[i]],value = T)),
        raster(grep("gamma_vh", dirs[[i]],value = T)),
        raster(grep("sigma_vv", dirs[[i]],value = T)),
        raster(grep("sigma_vh", dirs[[i]],value = T)))
  names(ind_stack) <- c("gamma_vv","gamma_vh","sigma_vv","sigma_vh")
  gamma_ratio <-ind_stack$gamma_vv/ind_stack$gamma_vh;names(gamma_ratio) <-"gamma_ratio"
  sigma_ratio <-ind_stack$sigma_vv/ind_stack$sigma_vh;names(sigma_ratio) <-"sigma_ratio"
  gamma_diff <-ind_stack$gamma_vv-ind_stack$gamma_vh;names(gamma_diff ) <- "gamma_diff"
  sigma_diff <-ind_stack$sigma_vv-ind_stack$sigma_vh;names(sigma_diff ) <- "sigma_diff"
  ind_stack <- stack(ind_stack,gamma_ratio, sigma_ratio, gamma_diff, sigma_diff)

  dat$gamma_ratio <- dat$gamma_vv/dat$gamma_vh
  dat$sigma_ratio <- dat$sigma_vv/dat$sigma_vh
  dat$gamma_diff <- dat$gamma_vv-dat$gamma_vh
  dat$sigma_diff <- dat$sigma_vv-dat$sigma_vh
  
  tmp <- glcm(raster(grep("gamma_vv", dirs[[i]],value = T)),na_opt="center",
              statistics = c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy", 
                  "second_moment"));names(tmp) <- gsub("glcm","gamma_vv",names(tmp))
  dat <- cbind(dat,raster::extract(tmp,sf::st_coordinates(trees)))
  ind_stack <- stack(ind_stack,tmp)
  tmp <- glcm(raster(grep("gamma_vh", dirs[[i]],value = T)),na_opt="center",
              statistics = c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy", 
                             "second_moment"));names(tmp) <- gsub("glcm","gamma_vh",names(tmp))
  dat <- cbind(dat,raster::extract(tmp,sf::st_coordinates(trees)))
  ind_stack <- stack(ind_stack,tmp)
  tmp <- glcm(raster(grep("sigma_vv", dirs[[i]],value = T)),na_opt="center",
              statistics = c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy", 
                             "second_moment"));names(tmp) <- gsub("glcm","sigma_vv",names(tmp))
  dat <- cbind(dat,raster::extract(tmp,sf::st_coordinates(trees)))
  ind_stack <- stack(ind_stack,tmp)
  tmp <- glcm(raster(grep("sigma_vh", dirs[[i]],value = T)),na_opt="center",
              statistics = c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy", 
                             "second_moment"));names(tmp) <- gsub("glcm","sigma_vh",names(tmp))
  dat <- cbind(dat,raster::extract(tmp,sf::st_coordinates(trees)))
  ind_stack <- stack(ind_stack,tmp)
  
  attributes(ind_stack)$meta <- meta
  list(dat,ind_stack)
},mc.cores = 1)

attributes(results[[10]][[1]])
sen1_res <- do.call("rbind", lapply(results,function(x) x[[1]]))
saveRDS(sen1_res, file="data/satellite_data/satellite_indices/sentinel1_indices.RDS")
saveRDS(results, file="data/satellite_data/satellite_results/prediction/sentinel1_predictors.RDS")

sen1_res_old <- readRDS("data/satellite_data/satellite_indices/sentinel1_indices_old.RDS")
results_old <- readRDS("data/satellite_data/satellite_results/prediction/sentinel1_predictors_old.RDS")
