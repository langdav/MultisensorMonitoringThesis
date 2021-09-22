#author: David Langenoh
#last modified: 17.08.2021
#description: script to extract NDVI time series per crown and fit logistic function to retrieve phenological start and end of budburst phase
#NOTE: no models fitted to Sentinel-2 data, as too little data available (3 days)

rm(list = ls())
library(rgdal);library(dplyr);library(RColorBrewer);library(ggplot2);library(ggpubr);library(lubridate);library(gridExtra)

#load data
load("out/log_function_models/all_values_fitted_models_output.RData")
load("out/log_function_models/all_values_fitted_models.RData")
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/mean_fitted_models.RData")
load("out/log_function_models/median_fitted_models_output.RData")
load("out/log_function_models/median_fitted_models.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models.RData")
# load("out/log_function_models/whole_forest_mean_fitted_models_output.RData")
# load("out/log_function_models/whole_forest_mean_fitted_models.RData")
# load("out/log_function_models/whole_forest_median_fitted_models_output.RData")
# load("out/log_function_models/whole_forest_median_fitted_models.RData")
# load("out/log_function_models/whole_forest_sentinel1_fitted_models_output.RData")
# load("out/log_function_models/whole_forest_sentinel1_fitted_models.RData")

budburst <- read.csv("data/budburst_data/budburst_long.csv")
budburst$date <- as.Date(budburst$date)

load("data/trees.RData")
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

#prepare data frames and functions
#create data frame for budburst percent vlines
buddies_budburst_percent <- NULL
for(tree in unique(budburst$tree_id)){
  tmp <- budburst %>% filter(tree_id == tree)
  percs <- unique(tmp$budburst_perc)
  for(perc in percs){
    if(perc == 0){
      buddies_budburst_percent <- rbind(buddies_budburst_percent, data.frame(tree_id = tree,
                                                                             first_date = yday(head(tmp$date,1)),
                                                                             last_date = yday(tmp$date[max(which(tmp$budburst_perc==0))]),
                                                                             budburst_perc = 0))
    }
    if(perc == 100){
      buddies_budburst_percent <- rbind(buddies_budburst_percent, data.frame(tree_id = tree,
                                                                             first_date = yday(tmp$date[min(which(tmp$budburst_perc==100))]),
                                                                             last_date = yday(tail(tmp$date,1)),
                                                                             budburst_perc = 100))
    }
    if(perc > 0 & perc < 100){
      buddies_budburst_percent <- rbind(buddies_budburst_percent, data.frame(tree_id = tree,
                                                                             first_date = yday(tmp$date[min(which(tmp$budburst_perc==perc))]),
                                                                             last_date = yday(tmp$date[max(which(tmp$budburst_perc==perc))]),
                                                                             budburst_perc = perc))
    }
  }
}

#create data frame for budburst start vline
buddies_budburst_start <- NULL
for(tree in unique(budburst$tree_id)){
  tmp <- budburst %>% filter(tree_id == tree)
  budburst_date_tmp <- yday(tmp$date[tmp$budburst == T][1])
  buddies_budburst_start <- rbind(buddies_budburst_start, data.frame(tree_id = tree,
                                                                     budburst_date = budburst_date_tmp))
}

#mean day of budburst from observations
budburst_mean_observation <- round(mean(buddies_budburst_start$budburst_date,na.rm=T),0)

#mean day of budburst from remote sensing
budburst_mean_remote <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  budburst_mean_mean <- round(mean(model_fitting_out_mean$SOS[model_fitting_out_mean$platform==platform],na.rm=T),0)
  budburst_mean_median <- round(mean(model_fitting_out_median$SOS[model_fitting_out_median$platform==platform],na.rm=T),0)
  budburst_mean_remote <- rbind(budburst_mean_remote, data.frame(platform = platform,
                                                                 mean_median = "mean",
                                                                 budburst_date = budburst_mean_mean,
                                                                 budburst_date_observation = budburst_mean_observation))
  budburst_mean_remote <- rbind(budburst_mean_remote, data.frame(platform = platform,
                                                                 mean_median = "median",
                                                                 budburst_date = budburst_mean_median,
                                                                 budburst_date_observation = budburst_mean_observation))
}
#add sentinel-1 derived budburst date
budburst_sen1_mean <- round(mean(model_fitting_out_sen1$SOS,na.rm=T),0)
budburst_mean_remote <- rbind(budburst_mean_remote, data.frame(platform = "sentinel1",
                                                               mean_median = NA,
                                                               budburst_date = budburst_sen1_mean,
                                                               budburst_date_observation = budburst_mean_observation))


#write.csv(budburst_mean_remote, "out/model_performance/mean_of_all_trees_budburst_dates.csv")

#plotting function; single tree and single platform
plot_SOS <- function(platform = "orthomosaic", mean_ndvi_values = T, median_ndvi_values = F, tree = 1, budburst_percent = F){
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
  
  #prepare data for plots
  
  # if(platform == "orthomosaic") {doylist <- c(113, 117, 134, 148, 155, 174)}
  # if(platform == "treetalker") {doylist <- c(96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 
  #                                            110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 
  #                                            123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 
  #                                            136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 
  #                                            149, 150, 151, 152)}
  # if(platform == "sentinel2") {doylist <- c(79, 89, 129, 164, 169)}
  # if(platform == "planetscope") {doylist <- c(66, 83, 89, 93, 99, 104, 110, 113, 114, 115, 116, 117, 128, 129, 144, 148, 151, 152)}
  # if(platform == "sentinel1") {doylist <- c(74, 77, 78, 80, 81, 83, 84, 86, 87, 89, 92, 93, 95, 96, 98, 99, 101, 102, 104, 105, 107, 
  #                                           108, 110, 111, 113, 114, 116, 117, 119, 122, 123, 125, 126, 128, 129, 131, 132, 134, 135, 
  #                                           137, 138, 140, 141, 143, 144, 146, 147, 149, 150)}
  
  
  ndvidat <- model_sel$ndvidat
  fitmodel <- model_sel$model
  if(!any(is.na(ndvidat)) & !any(is.na(fitmodel))){
    doylist <- model_sel$doylist
    platform <- model_sel$platform
    tree <- model_sel$tree
    preddoylist <- seq(head(doylist, 1),tail(doylist,1),1)
    predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
    
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
    
    
    if(budburst_percent == T){
      #plot with budburst percent
      buddies <- buddies_budburst_percent[buddies_budburst_percent$tree_id == tree_sel,]
      plot_out <- ggplot() +
        geom_point(aes(doylist, ndvidat)) +
        geom_line(aes(doylist,predNDVImod)) +
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
        geom_vline(xintercept = buddies$first_date[-1], linetype = "longdash", color = "red", size = .5) +
        geom_text(aes(x = buddies$first_date[-1],
                      y = max(ndvidat),
                      label = buddies$budburst_perc[-1],
                      vjust = 0,
                      angle = 90)) +
        xlim(xlim_min, xlim_max) +
        ylim(ylim_min, ylim_max) +
        theme_light() +
        xlab("Day of Year") +
        ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
        ggtitle(paste0(platform_sel,";",tree))
      
    } else {
      # plot with budburst start date
      buddies <- buddies_budburst_start[buddies_budburst_start$tree_id == tree_sel,]
      if(!is.na(data_sel$SOS)){
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_line(aes(doylist,predNDVImod)) +
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
          geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = buddies$budburst_date,
                        y = max(ndvidat),
                        label = "OBS",
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
          geom_line(aes(doylist,predNDVImod)) +
          geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = buddies$budburst_date,
                        y = max(ndvidat),
                        label = "OBS",
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
    
    
    if(budburst_percent == T){
      #plot with budburst percent
      buddies <- buddies_budburst_percent[buddies_budburst_percent$tree_id == tree_sel,]
      plot_out <- ggplot() +
        geom_point(aes(doylist, ndvidat)) +
        geom_vline(xintercept = buddies$first_date[-1], linetype = "longdash", color = "red", size = .5) +
        geom_text(aes(x = buddies$first_date[-1],
                      y = max(ndvidat),
                      label = buddies$budburst_perc[-1],
                      vjust = 0,
                      angle = 90)) +
        xlim(xlim_min, xlim_max) +
        ylim(ylim_min, ylim_max) +
        theme_light() +
        xlab("Day of Year") +
        ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
        ggtitle(paste0(platform_sel,";",tree))
      
    } else {
      # plot with budburst start date
      buddies <- buddies_budburst_start[buddies_budburst_start$tree_id == tree_sel,]
      if(!is.na(data_sel$SOS)){
        plot_out <- ggplot() +
          geom_point(aes(doylist, ndvidat)) +
          geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = buddies$budburst_date,
                        y = max(ndvidat),
                        label = "OBS",
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
          geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
          geom_text(aes(x = buddies$budburst_date,
                        y = max(ndvidat),
                        label = "OBS",
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
all_trees <- lapply(1:50, function(x) plot_SOS(platform = "sentinel2", mean_ndvi_values = F, median_ndvi_values = F, tree = x, budburst_percent = F))
all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
all_trees_panel


#Plotting: single tree and single platform
tree_no <- 1
one <- plot_SOS(platform = "orthomosaic", mean_ndvi_values = F, median_ndvi_values = F, tree = tree_no, budburst_percent = F)
two <- plot_SOS(platform = "planetscope", mean_ndvi_values = T, median_ndvi_values = F, tree = tree_no, budburst_percent = F)
three <- plot_SOS(platform = "treetalker", mean_ndvi_values = T, median_ndvi_values = F, tree = tree_no, budburst_percent = F)
four <- plot_SOS(platform = "sentinel2", mean_ndvi_values = T, median_ndvi_values = F, tree = tree_no, budburst_percent = F)
five <- plot_SOS(platform = "sentinel1", mean_ndvi_values = T, median_ndvi_values = F, tree = tree_no, budburst_percent = F)

ggpubr::ggarrange(one,two,three,four,five, ncol = 2, nrow = 3)



################################################################################

means_sen1 <- data.frame(platform = "sentinel1",
                         bb_mean = mean(model_fitting_out_sen1$SOS, na.rm = T)-124,
                         bb_sd = sd(model_fitting_out_sen1$SOS, na.rm = T),
                         ndvi = NA)

means_all <- NULL  
for(platform in unique(model_fitting_out_mean$platform)){
  tmp <- model_fitting_out_mean[which(model_fitting_out_mean$platform == platform),]
  m <- mean(tmp$SOS-124, na.rm = T)
  std <- sd(tmp$SOS, na.rm = T)
  
  means_all <- rbind(means_all, data.frame(platform = platform,
                                           bb_mean = m,
                                           bb_sd = std,
                                           ndvi = "mean"))
}

for(platform in unique(model_fitting_out_median$platform)){
  tmp <- model_fitting_out_median[which(model_fitting_out_median$platform == platform),]
  m <- mean(tmp$SOS-124, na.rm = T)
  std <- sd(tmp$SOS, na.rm = T)
  
  means_all <- rbind(means_all, data.frame(platform = platform,
                                           bb_mean = m,
                                           bb_sd = std,
                                           ndvi = "median"))
}

for(platform in unique(model_fitting_out_all$platform)){
  tmp <- model_fitting_out_all[which(model_fitting_out_all$platform == platform),]
  m <- mean(tmp$SOS-124, na.rm = T)
  std <- sd(tmp$SOS, na.rm = T)
  
  means_all <- rbind(means_all, data.frame(platform = platform,
                                           bb_mean = m,
                                           bb_sd = std,
                                           ndvi = "all_values"))
}
means_all <- rbind(means_all, means_sen1)



ggplot(data = means_all, aes(x = bb_mean, y = platform, colour = ndvi)) + 
  geom_point() +
  geom_errorbarh(aes(xmin=bb_mean-bb_sd, xmax=bb_mean+bb_sd)) +
  geom_vline(xintercept = 0, linetype = "longdash", color = "red", size = .5) +
  facet_grid()

























# #plotting function; all trees and single platform
# plot_SOS_whole_forest <- function(platform = "orthomosaic", mean_ndvi_values = T){
#   if(platform == "sentinel1"){
#     model_fitting_out <- whole_forest_model_fitting_out_sen1
#     models <- whole_forest_models_sen1
#   } else {
#     if(mean_ndvi_values == T){
#       model_fitting_out <- whole_forest_model_fitting_out_mean
#       models <- whole_forest_models_mean
#     } else {
#       model_fitting_out <- whole_forest_model_fitting_out_median
#       models <- whole_forest_models_median
#     } 
#   }
#   
#   #select platform; one of sentinel2, orthomosaic, treetalker, planetscope
#   platform_sel <- platform
#   
#   #grab a model and respective SOS, MOS, EOS data
#   if(platform == "sentinel1"){
#     model_sel <- models[[1]]
#   } else {
#     entry <- NULL
#     for(i in 1:length(models)){
#       if(models[[i]][["platform"]] == platform_sel){
#         entry <- c(entry, i)
#       }
#     }
#     
#     model_sel <- models[[entry]]
#   }
#   
#   if(platform == "sentinel1"){
#     data_sel <- model_fitting_out
#   } else {
#     data_sel <- model_fitting_out[which(model_fitting_out$platform == platform_sel),]
#   }  
#   
#   #prepare data for plots
#   ndvidat <- model_sel$ndvidat
#   doylist <- model_sel$doylist
#   fitmodel <- model_sel$model
#   platform <- model_sel$platform
#   #tree <- model_sel$tree
#   preddoylist <- seq(head(doylist, 1),tail(doylist,1),1)
#   predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
#   
#   #x- and y-axis limits
#   if(platform == "sentinel1"){
#     load("out/all_in_one/sentinel1_daily_sigma_ratio_all_trees.RData")
#     xlim_min <- min(unique(sen1_all_tree_means$doy))
#     xlim_max <- max(unique(sen1_all_tree_means$doy))
#     ylim_min <- min(unique(sen1_all_tree_means$sigma_ratio_mean))
#     ylim_max <- max(unique(sen1_all_tree_means$sigma_ratio_mean))
#   } else {
#     load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
#     aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
#     aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
#     
#     xlim_min <- yday(min(unique(aio$date)))
#     xlim_max <- yday(max(unique(aio$date)))
#     ylim_min <- min(unique(aio$ndvi_mean))
#     ylim_max <- max(unique(aio$ndvi_mean))
#   }
#   
#   # plot with budburst start date
#   
#   plot_out <- ggplot() +
#     geom_point(aes(doylist, ndvidat)) +
#     geom_line(aes(doylist,predNDVImod)) +
#     #geom_vline(xintercept = c(data_sel$SOS,data_sel$MOS,data_sel$EOS), linetype = "dotdash", color = "blue", size = .8) +
#     geom_vline(xintercept = data_sel$SOS, linetype = "dotdash", color = "blue", size = .8) +
#     # geom_text(aes(x = c(data_sel$SOS,data_sel$MOS,data_sel$EOS),
#     #               y = min(ndvidat),
#     #               label = c("SOS","MOS","EOS"),
#     #               vjust = 0,
#     #               angle = 90)) +
#     geom_text(aes(x = data_sel$SOS,
#                   y = min(ndvidat),
#                   label = "predicted",
#                   vjust = 0,
#                   angle = 90)) +
#     geom_vline(xintercept = budburst_mean_observation, linetype = "longdash", color = "red", size = .5) +
#     geom_text(aes(x = budburst_mean_observation,
#                   y = max(ndvidat),
#                   label = "OBS",
#                   vjust = 0,
#                   angle = 90)) +
#     xlim(xlim_min, xlim_max) +
#     ylim(ylim_min, ylim_max) +
#     theme_light() +
#     xlab("Day of Year") +
#     ylab(ifelse(platform == "sentinel1","Backscatter/db","NDVI")) +
#     ggtitle(platform_sel)
#   
#   return(plot_out)
# }
# 
# #plot mean budburst date over all trees
# plot_SOS_whole_forest(platform = "sentinel1", mean_ndvi_values = T)
# 
# one <- plot_SOS_whole_forest(platform = "orthomosaic", mean_ndvi_values = T)
# two <- plot_SOS_whole_forest(platform = "planetscope", mean_ndvi_values = T)
# three <- plot_SOS_whole_forest(platform = "treetalker", mean_ndvi_values = T)
# four <- plot_SOS_whole_forest(platform = "sentinel2", mean_ndvi_values = T)
# five  <- plot_SOS_whole_forest(platform = "sentinel1", mean_ndvi_values = T)
# ggpubr::ggarrange(one,two,three,four,five, ncol = 2, nrow = 3)
# 
# one <- plot_SOS_whole_forest(platform = "orthomosaic", mean_ndvi_values = F)
# two <- plot_SOS_whole_forest(platform = "planetscope", mean_ndvi_values = F)
# three <- plot_SOS_whole_forest(platform = "treetalker", mean_ndvi_values = F)
# four <- plot_SOS_whole_forest(platform = "sentinel2", mean_ndvi_values = F)
# five  <- plot_SOS_whole_forest(platform = "sentinel1", mean_ndvi_values = F)
# ggpubr::ggarrange(one,two,three,four,five, ncol = 2, nrow = 3)
