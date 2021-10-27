#author: David Langenoh
#last modified: 17.08.2021
#description: script to extract NDVI time series per crown and fit logistic function to retrieve phenological start and end of budburst phase
#NOTE: no models fitted to Sentinel-2 data, as too little data available (3 days)

rm(list = ls())
library(rgdal);library(dplyr);library(RColorBrewer);library(ggplot2);library(ggpubr);library(lubridate);library(gridExtra)

#load data
# load("out/log_function_models/all_values_fitted_models_output.RData")
# load("out/log_function_models/all_values_fitted_models.RData")
# load("out/log_function_models/median_fitted_models_output.RData")
# load("out/log_function_models/median_fitted_models.RData")
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/mean_fitted_models.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models.RData")

load("data/trees.RData")
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

#prepare data frames and functions
#create data frame for budburst start vline
# budburst <- read.csv("data/budburst_data/budburst_long.csv")
# budburst$date <- as.Date(budburst$date)
# 
# buddies_budburst_start <- NULL
# for(tree in unique(budburst$tree_id)){
#   tmp <- budburst %>% filter(tree_id == tree)
#   budburst_date_tmp <- yday(tmp$date[tmp$budburst == T][1])
#   buddies_budburst_start <- rbind(buddies_budburst_start, data.frame(tree_id = tree,
#                                                                      budburst_date = budburst_date_tmp))
# }

#plotting function; single tree and single platform
plot_SOS <- function(platform = "orthomosaic", obs_phase = "phase_d", mean_ndvi_values = T, median_ndvi_values = F, tree = 1, budburst_percent = F){
  if(platform == "sentinel1"){
    models <- models_sen1
    model_fitting_out <- model_fitting_out_sen1
  } else {
    if(mean_ndvi_values == T & median_ndvi_values == F){
      models <- models_mean
      model_fitting_out <- model_fitting_out_mean
    } else if(mean_ndvi_values == F & median_ndvi_values == T){
      models <- models_median
      model_fitting_out <- model_fitting_out_median
    } else {
      models <- models_all
      model_fitting_out <- model_fitting_out_all
    }
  }
  
  #select platform; one of sentinel2, orthomosaic, treetalker, planetscope
  platform_sel <- platform
  
  #select tree; between 1 and 50
  tree_sel <- trees$tree_id[tree]
  
  #grab a model and respective SOS, MOS, EOS data
  if(platform == "sentinel1"){
    entry <- NULL
    for(i in 1:length(models)){
      if(models[[i]][["tree_id"]] == tree_sel){
        entry <- c(entry, i)
      }
    }
    model_sel <- models[[entry]]
  } else {
    entry <- NULL
    for(i in 1:length(models)){
      if(models[[i]][["tree_id"]] == tree_sel & models[[i]][["platform"]] == platform_sel){
        entry <- c(entry, i)
      }
    }
    
    model_sel <- models[[entry]]
  }
  
  if(platform == "sentinel1"){
    data_sel <- model_fitting_out[which(model_fitting_out$tree_id == tree_sel),]
  } else {
    data_sel <- model_fitting_out[which(model_fitting_out$platform == platform_sel & model_fitting_out$tree_id == tree_sel),]
  }  
  
  ndvidat <- model_sel$ndvidat
  ndvi_sd <- model_sel[[4]]
  fitmodel <- model_sel$model
  if(!any(is.na(ndvidat)) & !any(is.na(fitmodel))){
    doylist <- model_sel$doylist
    platform <- model_sel$platform
    tree <- model_sel$tree
    preddoylist <- seq(head(doylist, 1),tail(doylist,1),1)
    predNDVImod <- predict(fitmodel,data.frame(doylist=preddoylist))
    
    #x- and y-axis limits
    if(platform == "sentinel1"){
      load("out/sentinel1/backscatter_all_with_phenoclasses_sentinel1.RData")
      xlim_min <- min(unique(sen1_backscatter$doy))
      xlim_max <- max(unique(sen1_backscatter$doy))
      ylim_min <- min(unique(sen1_backscatter$sigma_ratio))
      ylim_max <- max(unique(sen1_backscatter$sigma_ratio))
    } else {
      if(mean_ndvi_values == T & median_ndvi_values == F){
        load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
        aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
        aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
        
        xlim_min <- yday(min(unique(aio$date)))
        xlim_max <- yday(max(unique(aio$date)))
        ylim_min <- min(unique(aio$ndvi_mean))
        ylim_max <- max(unique(aio$ndvi_mean))
      } else if(mean_ndvi_values == F & median_ndvi_values == T){
        load("out/all_in_one/aio_daily_ndvi_per_tree_medians.RData")
        aio <- aio_daily_ndvi_medians;rm(aio_daily_ndvi_medians)
        aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
        
        xlim_min <- yday(min(unique(aio$date)))
        xlim_max <- yday(max(unique(aio$date)))
        ylim_min <- min(unique(aio$ndvi_median))
        ylim_max <- max(unique(aio$ndvi_median))
      } else {
        load("out/all_in_one/aio_daily_ndvi_all.RData")
        aio <- aio_daily_ndvi_all;rm(aio_daily_ndvi_all)
        aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
        
        xlim_min <- yday(min(unique(aio$date)))
        xlim_max <- yday(max(unique(aio$date)))
        ylim_min <- min(unique(aio$ndvi))
        ylim_max <- max(unique(aio$ndvi))
      }
    }
    
    # buddies <- buddies_budburst_start[buddies_budburst_start$tree_id == tree_sel,]
    if(mean_ndvi_values == T & median_ndvi_values == F){
      if(!is.na(data_sel$SOS)){
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_errorbar(aes(x = doylist ,ymin=ndvidat-ndvi_sd, ymax=ndvidat+ndvi_sd), width = 1) +
          geom_line(aes(preddoylist,predNDVImod)) +
          #geom_vline(xintercept = c(data_sel$SOS,data_sel$MOS,data_sel$EOS), linetype = "dotdash", color = "blue", size = .8) +
          geom_vline(xintercept = data_sel$SOS, linetype = "dotdash", color = "blue", size = .8) +
          # geom_text(aes(x = c(data_sel$SOS,data_sel$MOS,data_sel$EOS),
          #               y = min(ndvidat),
          #               label = c("SOS","MOS","EOS"),
          #               vjust = 0,
          #               angle = 90)) +
          geom_text(aes(x = data_sel$SOS,
                        y = min(ndvidat),
                        label = "predicted",
                        vjust = 0,
                        angle = 90)) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          # xlim(xlim_min, xlim_max) +
          # ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      } else {
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_errorbar(aes(x = doylist ,ymin=ndvidat-ndvi_sd, ymax=ndvidat+ndvi_sd), width = 1) +
          geom_line(aes(preddoylist,predNDVImod)) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          xlim(xlim_min, xlim_max) +
          ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      }
    } else {
      if(!is.na(data_sel$SOS)){
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_line(aes(preddoylist,predNDVImod)) +
          #geom_vline(xintercept = c(data_sel$SOS,data_sel$MOS,data_sel$EOS), linetype = "dotdash", color = "blue", size = .8) +
          geom_vline(xintercept = data_sel$SOS, linetype = "dotdash", color = "blue", size = .8) +
          # geom_text(aes(x = c(data_sel$SOS,data_sel$MOS,data_sel$EOS),
          #               y = min(ndvidat),
          #               label = c("SOS","MOS","EOS"),
          #               vjust = 0,
          #               angle = 90)) +
          geom_text(aes(x = data_sel$SOS,
                        y = min(ndvidat),
                        label = "predicted",
                        vjust = 0,
                        angle = 90)) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          xlim(xlim_min, xlim_max) +
          ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      } else {
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_line(aes(preddoylist,predNDVImod)) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          xlim(xlim_min, xlim_max) +
          ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      }
    }
    
    
    
  } else if(!any(is.na(ndvidat)) & any(is.na(fitmodel))){
    doylist <- model_sel$doylist
    platform <- model_sel$platform
    tree <- model_sel$tree
    
    #x- and y-axis limits
    if(platform == "sentinel1"){
      load("out/sentinel1/backscatter_all_with_phenoclasses_sentinel1.RData")
      xlim_min <- min(unique(sen1_backscatter$doy))
      xlim_max <- max(unique(sen1_backscatter$doy))
      ylim_min <- min(unique(sen1_backscatter$sigma_ratio))
      ylim_max <- max(unique(sen1_backscatter$sigma_ratio))
    } else {
      if(mean_ndvi_values == T & median_ndvi_values == F){
        load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
        aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
        aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
        
        xlim_min <- yday(min(unique(aio$date)))
        xlim_max <- yday(max(unique(aio$date)))
        ylim_min <- min(unique(aio$ndvi_mean))
        ylim_max <- max(unique(aio$ndvi_mean))
      } else if(mean_ndvi_values == F & median_ndvi_values == T){
        load("out/all_in_one/aio_daily_ndvi_per_tree_medians.RData")
        aio <- aio_daily_ndvi_medians;rm(aio_daily_ndvi_medians)
        aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
        
        xlim_min <- yday(min(unique(aio$date)))
        xlim_max <- yday(max(unique(aio$date)))
        ylim_min <- min(unique(aio$ndvi_median))
        ylim_max <- max(unique(aio$ndvi_median))
      } else {
        load("out/all_in_one/aio_daily_ndvi_all.RData")
        aio <- aio_daily_ndvi_all;rm(aio_daily_ndvi_all)
        aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
        
        xlim_min <- yday(min(unique(aio$date)))
        xlim_max <- yday(max(unique(aio$date)))
        ylim_min <- min(unique(aio$ndvi))
        ylim_max <- max(unique(aio$ndvi))
      }
    }
    
    # buddies <- buddies_budburst_start[buddies_budburst_start$tree_id == tree_sel,]
    if(mean_ndvi_values == T & median_ndvi_values == F){
      if(!is.na(data_sel$SOS)){
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_errorbar(aes(x = doylist ,ymin=ndvidat-ndvi_sd, ymax=ndvidat+ndvi_sd), width = 1) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          xlim(xlim_min, xlim_max) +
          ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      } else {
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_errorbar(aes(x = doylist ,ymin=ndvidat-ndvi_sd, ymax=ndvidat+ndvi_sd), width = 1) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          xlim(xlim_min, xlim_max) +
          ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      }
    } else {
      if(!is.na(data_sel$SOS)){
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          xlim(xlim_min, xlim_max) +
          ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      } else {
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          # geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = buddies$budburst_date,
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          # geom_vline(xintercept = data_sel[,paste0("SOS_",obs_phase)], linetype = "longdash", color = "red", size = .5) +
          # geom_text(aes(x = data_sel[,paste0("SOS_",obs_phase)],
          #               y = max(ndvidat),
          #               label = "OBS",
          #               vjust = 0,
          #               angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_d, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_d,
                        y = max(ndvidat),
                        label = "phase_d",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_e, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_e,
                        y = max(ndvidat),
                        label = "phase_e",
                        vjust = 0,
                        angle = 90)) +
          geom_vline(xintercept = data_sel$SOS_phase_f, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = data_sel$SOS_phase_f,
                        y = max(ndvidat),
                        label = "phase_f",
                        vjust = 0,
                        angle = 90)) +
          xlim(xlim_min, xlim_max) +
          ylim(ylim_min, ylim_max) +
          theme_light() +
          xlab("Day of Year") +
          ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
          ggtitle(paste0(platform_sel,";",tree))
      }
    }
    
  } else {
    plot_out <- ggplot() +
      geom_text(aes(x = 1, y = 1, label = "no data available"), size = 5) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.background = element_rect(fill = 'white', colour = 'white')) +
      ggtitle(paste0(platform_sel,";",tree_sel))
  }
  
  return(plot_out)
}

#Plotting: all trees, single platform; one panel
# for(platform in c(unique(model_fitting_out_all$platform),"sentinel1")){
#   for(phase in c("phase_d","phase_e","phase_f")){
#     cat("Processing", platform, ",", phase,"\n")
#     
#     all_trees <- lapply(1:25, function(x) plot_SOS(platform = platform, obs_phase = phase, mean_ndvi_values = T, median_ndvi_values = F, tree = x, budburst_percent = F))
#     all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
#     ggsave(paste0("out/model_plots/",platform,"_mean_",phase,"_1_25",".png"),all_trees_panel,width = 20,height=12)
#     
#     all_trees <- lapply(26:50, function(x) plot_SOS(platform = platform, obs_phase = phase, mean_ndvi_values = T, median_ndvi_values = F, tree = x, budburst_percent = F))
#     all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
#     ggsave(paste0("out/model_plots/",platform,"_mean_",phase,"_25_50",".png"),all_trees_panel,width = 20,height=12)
#     
#     all_trees <- lapply(1:25, function(x) plot_SOS(platform = platform, obs_phase = phase, mean_ndvi_values = F, median_ndvi_values = T, tree = x, budburst_percent = F))
#     all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
#     ggsave(paste0("out/model_plots/",platform,"_median_",phase,"_1_25",".png"),all_trees_panel,width = 20,height=12)
#     
#     all_trees <- lapply(26:50, function(x) plot_SOS(platform = platform, obs_phase = phase, mean_ndvi_values = F, median_ndvi_values = T, tree = x, budburst_percent = F))
#     all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
#     ggsave(paste0("out/model_plots/",platform,"_median_",phase,"_25_50",".png"),all_trees_panel,width = 20,height=12)
#     
#     all_trees <- lapply(1:25, function(x) plot_SOS(platform = platform, obs_phase = phase, mean_ndvi_values = F, median_ndvi_values = F, tree = x, budburst_percent = F))
#     all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
#     ggsave(paste0("out/model_plots/",platform,"_all_values_",phase,"_1_25",".png"),all_trees_panel,width = 20,height=12)
#     
#     all_trees <- lapply(26:50, function(x) plot_SOS(platform = platform, obs_phase = phase, mean_ndvi_values = F, median_ndvi_values = F, tree = x, budburst_percent = F))
#     all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
#     ggsave(paste0("out/model_plots/",platform,"_all_values_",phase,"_25_50",".png"),all_trees_panel,width = 20,height=12)
#   }
# }
for(platform in c(unique(model_fitting_out_mean$platform),"sentinel1")){
    cat("Processing", platform,"\n")
    
    all_trees <- lapply(1:25, function(x) plot_SOS(platform = platform, obs_phase = "phase_d", mean_ndvi_values = T, median_ndvi_values = F, tree = x, budburst_percent = F))
    all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
    ggsave(paste0("out/model_plots/",platform,"_mean_1_25",".png"),all_trees_panel,width = 20,height=12)
    
    all_trees <- lapply(26:50, function(x) plot_SOS(platform = platform, obs_phase = "phase_d", mean_ndvi_values = T, median_ndvi_values = F, tree = x, budburst_percent = F))
    all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
    ggsave(paste0("out/model_plots/",platform,"_mean_25_50",".png"),all_trees_panel,width = 20,height=12)
    
    # all_trees <- lapply(1:25, function(x) plot_SOS(platform = platform, obs_phase = "phase_d", mean_ndvi_values = F, median_ndvi_values = T, tree = x, budburst_percent = F))
    # all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
    # ggsave(paste0("out/model_plots/",platform,"_median_1_25",".png"),all_trees_panel,width = 20,height=12)
    # 
    # all_trees <- lapply(26:50, function(x) plot_SOS(platform = platform, obs_phase = "phase_d", mean_ndvi_values = F, median_ndvi_values = T, tree = x, budburst_percent = F))
    # all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
    # ggsave(paste0("out/model_plots/",platform,"_median_25_50",".png"),all_trees_panel,width = 20,height=12)
    # 
    # all_trees <- lapply(1:25, function(x) plot_SOS(platform = platform, obs_phase = "phase_d", mean_ndvi_values = F, median_ndvi_values = F, tree = x, budburst_percent = F))
    # all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
    # ggsave(paste0("out/model_plots/",platform,"_all_values_1_25",".png"),all_trees_panel,width = 20,height=12)
    # 
    # all_trees <- lapply(26:50, function(x) plot_SOS(platform = platform, obs_phase = "phase_d", mean_ndvi_values = F, median_ndvi_values = F, tree = x, budburst_percent = F))
    # all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
    # ggsave(paste0("out/model_plots/",platform,"_all_values_25_50",".png"),all_trees_panel,width = 20,height=12)
}

all_trees <- lapply(1:50, function(x) plot_SOS(platform = "planetscope", mean_ndvi_values = T, median_ndvi_values = F, tree = x, budburst_percent = F))
all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
all_trees_panel
#ggsave("out/model_plots/test.png",all_trees_panel)

#Plotting: single tree and single platform
tree_no <- 3
plot_SOS(platform = "sentinel1", obs_phase = "phase_d", mean_ndvi_values = T, median_ndvi_values = F, tree = tree_no, budburst_percent = F)

#create plot of each platforms deviation from mean budburst date (with errorbars)
# mean observed BB date
load("out/log_function_models/budburst_fitted_models_output.RData")
mean_obs_d <- mean(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == "phase_d")], na.rm = T)
sd_obs_d <- sd(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == "phase_d")], na.rm = T)

mean_obs_e <- mean(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == "phase_e")], na.rm = T)
sd_obs_e <- sd(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == "phase_e")], na.rm = T)

mean_obs_f <- mean(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == "phase_f")], na.rm = T)
sd_obs_f <- sd(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == "phase_f")], na.rm = T)


# mean predicted dates per platform
means_all <- NULL  
for(platform in unique(model_fitting_out_mean$platform)){
  tmp <- model_fitting_out_mean[which(model_fitting_out_mean$platform == platform),]
  m_d <- mean(tmp$SOS-mean_obs_d, na.rm = T)
  m_e <- mean(tmp$SOS-mean_obs_e, na.rm = T)
  m_f <- mean(tmp$SOS-mean_obs_f, na.rm = T)
  std <- sd(tmp$SOS, na.rm = T)
  
  means_all <- rbind(means_all, data.frame(platform = platform,
                                           bb_mean_d = m_d,
                                           bb_mean_e = m_e,
                                           bb_mean_f = m_f,
                                           bb_sd = std,
                                           ndvi = "mean"))
  
  tmp <- model_fitting_out_median[which(model_fitting_out_median$platform == platform),]
  m_d <- mean(tmp$SOS-mean_obs_d, na.rm = T)
  m_e <- mean(tmp$SOS-mean_obs_e, na.rm = T)
  m_f <- mean(tmp$SOS-mean_obs_f, na.rm = T)
  std <- sd(tmp$SOS, na.rm = T)
  
  means_all <- rbind(means_all, data.frame(platform = platform,
                                           bb_mean_d = m_d,
                                           bb_mean_e = m_e,
                                           bb_mean_f = m_f,
                                           bb_sd = std,
                                           ndvi = "median"))
  
  tmp <- model_fitting_out_all[which(model_fitting_out_all$platform == platform),]
  m_d <- mean(tmp$SOS-mean_obs_d, na.rm = T)
  m_e <- mean(tmp$SOS-mean_obs_e, na.rm = T)
  m_f <- mean(tmp$SOS-mean_obs_f, na.rm = T)
  std <- sd(tmp$SOS, na.rm = T)
  
  means_all <- rbind(means_all, data.frame(platform = platform,
                                           bb_mean_d = m_d,
                                           bb_mean_e = m_e,
                                           bb_mean_f = m_f,
                                           bb_sd = std,
                                           ndvi = "all_values"))
}
means_sen1 <- data.frame(platform = "sentinel1",
                         bb_mean_d = mean(model_fitting_out_sen1$SOS, na.rm = T)-mean_obs_d,
                         bb_mean_e = mean(model_fitting_out_sen1$SOS, na.rm = T)-mean_obs_e,
                         bb_mean_f = mean(model_fitting_out_sen1$SOS, na.rm = T)-mean_obs_f,
                         bb_sd = sd(model_fitting_out_sen1$SOS, na.rm = T),
                         ndvi = "mean")

load("out/log_function_models/multiplatform_fitted_models_output.RData")
means_multi <- data.frame(platform = "multiplatform",
                         bb_mean_d = mean(model_fitting_out_multiplatform$SOS, na.rm = T)-mean_obs_d,
                         bb_mean_e = mean(model_fitting_out_multiplatform$SOS, na.rm = T)-mean_obs_e,
                         bb_mean_f = mean(model_fitting_out_multiplatform$SOS, na.rm = T)-mean_obs_f,
                         bb_sd = sd(model_fitting_out_multiplatform$SOS, na.rm = T),
                         ndvi = "mean")

means_all <- rbind(means_all, means_sen1,means_multi)


means_all %>% 
  filter(ndvi == "mean") %>% 
  ggplot(aes(x = bb_mean_d, y = platform, color = platform)) + 
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin=bb_mean_d-bb_sd, xmax=bb_mean_d+bb_sd), height = .3) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "black", size = .5) +
  geom_vline(xintercept = sd_obs_d, linetype = "dashed", color = "black", size = .5) +
  geom_vline(xintercept = 0-sd_obs_d, linetype = "dashed", color = "black", size = .5) +
  #facet_grid(ndvi ~ .) +
  theme_light() +
  xlab("Deviation from observed prediction (Phase D)") +
  ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
  ggtitle("per-platform deviation (in days) from predicted budbust date from observed values (Phase D)")
ggsave("out/model_plots/method_deviations_from_obs_d.png",width = 10,height=6)

means_all %>% 
  filter(ndvi == "mean") %>% 
  ggplot(aes(x = bb_mean_e, y = platform, color = platform)) + 
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin=bb_mean_e-bb_sd, xmax=bb_mean_e+bb_sd), height = .3) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "black", size = .5) +
  geom_vline(xintercept = sd_obs_e, linetype = "dashed", color = "black", size = .5) +
  geom_vline(xintercept = 0-sd_obs_e, linetype = "dashed", color = "black", size = .5) +
  #facet_grid(ndvi ~ .) +
  theme_light() +
  xlab("Deviation from observed prediction (Phase E)") +
  ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
  ggtitle("per-platform deviation (in days) from predicted budbust date from observed values (Phase E)")
ggsave("out/model_plots/method_deviations_from_obs_e.png",width = 10,height=6)

means_all %>% 
  filter(ndvi == "mean") %>% 
  ggplot(aes(x = bb_mean_f, y = platform, color = platform)) + 
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin=bb_mean_f-bb_sd, xmax=bb_mean_f+bb_sd), height = .3) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "black", size = .5) +
  geom_vline(xintercept = sd_obs_f, linetype = "dashed", color = "black", size = .5) +
  geom_vline(xintercept = 0-sd_obs_f, linetype = "dashed", color = "black", size = .5) +
  #facet_grid(ndvi ~ .) +
  theme_light() +
  xlab("Deviation from observed prediction (Phase F)") +
  ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
  ggtitle("per-platform deviation (in days) from predicted budbust date from observed values (Phase F)")
ggsave("out/model_plots/method_deviations_from_obs_f.png",width = 10,height=6)

######################################################
means_all <- means_all[means_all$ndvi == "mean",]
test <- data.frame(platform = rep(means_all$platform,3),
                   mean = c(means_all$bb_mean_d,means_all$bb_mean_e,means_all$bb_mean_f),
                   sd = rep(means_all$bb_sd,3),
                   phase = c(rep("phase_d",5),rep("phase_e",5),rep("phase_f",5)))
test %>% 
  ggplot(aes(x = mean, y = platform, color = platform)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin=mean-sd, xmax=mean+sd), height = .3) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "black", size = .5) +
  geom_vline(xintercept = sd_obs_f, linetype = "dashed", color = "black", size = .5) +
  geom_vline(xintercept = 0-sd_obs_f, linetype = "dashed", color = "black", size = .5) +
  facet_grid(phase ~ .) +
  theme_light() +
  xlab("Deviation: Predicted - Observed") +
  ylab("Platform")
ggsave("out/model_plots/method_deviations_from_obs_mean_only.png",width = 10,height=6)
