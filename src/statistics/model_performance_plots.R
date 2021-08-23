#author: David Langenohl; based on the work of Dominic Fawcett
#last modified: 17.08.2021
#description: script to extract NDVI time series per crown and fit logistic function to retrieve phenological start and end of budburst phase
#NOTE: no models fitted to Sentinel-2 data, as too little data available (3 days)

rm(list = ls())
library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);library(minpack.lm) #minpack.lm isneeded for the use of minpack.lm::nlsLM()

#load data
load("out/log_function_models/fitted_models_output.RData")
load("out/log_function_models/fitted_models.RData")

budburst <- read.csv("data/budburst_data/budburst_long.csv")
budburst$date <- as.Date(budburst$date)

load("data/trees.RData")
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

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
buddies_budburs_start <- NULL
for(tree in unique(budburst$tree_id)){
  tmp <- budburst %>% filter(tree_id == tree)
  budburst_date_tmp <- yday(tmp$date[tmp$budburst == T][1])
  buddies_budburs_start <- rbind(buddies_budburs_start, data.frame(tree_id = tree,
                                                                   budburst_date = budburst_date_tmp))
}

#plotting function
plot_SOS <- function(platform = "orthomosaic", tree = 1, budburst_percent = F){
  #select platform; one of sentinel2, orthomosaic, treetalker, planetscope
  platform_sel <- platform
  
  #select tree; between 1 and 50
  tree_sel <- trees$tree_id[tree]
  
  #grab a model and respective SOS, MOS, EOS data
  entry <- NULL
  for(i in 1:length(models)){
    if(models[[i]][["tree"]] == tree_sel & models[[i]][["platform"]] == platform_sel){
      entry <- c(entry, i)
    }
  }
  
  model_sel <- models[[entry]]
  data_sel <- model_fitting_out[which(model_fitting_out$platform == platform_sel & model_fitting_out$tree == tree_sel),]
  
  #prepare data for plots
  ndvidat <- model_sel$ndvidat
  doylist <- model_sel$doylist
  fitmodel <- model_sel$model
  platform <- model_sel$platform
  tree <- model_sel$tree
  preddoylist <- seq(head(doylist, 1),tail(doylist,1),1)
  predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
  
  #x- and y-axis limits
  load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
  aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
  aio <- aio[-which(aio$date > as.Date("2021-07-01")),]
  
  xlim_min <- yday(min(unique(aio$date)))
  xlim_max <- yday(max(unique(aio$date)))
  ylim_min <- min(unique(aio$ndvi_mean))
  ylim_max <- max(unique(aio$ndvi_mean))
  
  
  if(budburst_percent == T){
    #plot with budburst percent
    buddies <- buddies_budburst_percent[buddies_budburst_percent$tree_id == tree_sel,]
    plot_out <- ggplot() +
      geom_point(aes(doylist, ndvidat)) +
      geom_line(aes(doylist,predNDVImod)) +
      geom_vline(xintercept = c(data_sel$SOS,data_sel$MOS,data_sel$EOS), linetype = "dotdash", color = "blue", size = .8) +
      geom_text(aes(x = c(data_sel$SOS,data_sel$MOS,data_sel$EOS),
                    y = min(ndvidat),
                    label = c("SOS","MOS","EOS"),
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
      ylab("NDVI")
  } else {
    # plot with budburst start date
    buddies <- buddies_budburs_start[buddies_budburs_start$tree_id == tree_sel,]
    plot_out <- ggplot() +
      geom_point(aes(doylist, ndvidat)) +
      geom_line(aes(doylist,predNDVImod)) +
      geom_vline(xintercept = c(data_sel$SOS,data_sel$MOS,data_sel$EOS), linetype = "dotdash", color = "blue", size = .8) +
      geom_text(aes(x = c(data_sel$SOS,data_sel$MOS,data_sel$EOS),
                    y = min(ndvidat),
                    label = c("SOS","MOS","EOS"),
                    vjust = 0,
                    angle = 90)) +
      geom_vline(xintercept = buddies$budburst_date, linetype = "longdash", color = "red", size = .5) +
      geom_text(aes(x = buddies$budburst_date,
                    y = max(ndvidat),
                    label = "budburst",
                    vjust = 0,
                    angle = 90)) +
      xlim(xlim_min, xlim_max) +
      ylim(ylim_min, ylim_max) +
      theme_light() +
      xlab("Day of Year") +
      ylab("NDVI")
  }
  return(plot_out)
}

###############################################
#Plotting: single tree and single platform
# budbust start only
plot_SOS(platform = "orthomosaic", tree = 1, budburst_percent = F)

# budbust percent
plot_SOS(platform = "orthomosaic", tree = 1, budburst_percent = T)

######
one <- plot_SOS(platform = "orthomosaic", tree = 1, budburst_percent = F)
two <- plot_SOS(platform = "planetscope", tree = 1, budburst_percent = F)
three <- plot_SOS(platform = "treetalker", tree = 1, budburst_percent = F)
four <- plot_SOS(platform = "sentinel2", tree = 1, budburst_percent = F)

library(ggpubr)
ggpubr::ggarrange(one,two,three,four, ncol = 2, nrow = 2)





