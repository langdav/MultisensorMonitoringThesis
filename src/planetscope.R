library(dplyr);library(ggplot2);library(viridis);library(tidyr)
library(raster);library(rgeos);library(rgdal);library(lidR)

## load data
# phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

# trees
trees <- readRDS("data/trees.RDS")

## create colorblind friendle palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#############################################################
## create NDVI (or other indices) from planetscope images ##
###########################################################
## NDVI for all planetscope files; create raster stack
files_planetscope <- list.files("data/satellite_data/planetscope/",pattern = "SR_clip.tif")
for(i in files_planetscope){
  tmp_stack <- stack(paste0("data/satellite_data/planetscope/", i)) #import data
  #ndvi <- (tmp_stack[[4]] - tmp_stack[[3]])/(tmp_stack[[4]] + tmp_stack[[3]]) NDVI = (nir-red)/(nir+red)
  ndvi <- RStoolbox::spectralIndices(tmp_stack, blue = 1, green = 2, red = 3, nir = 4, indices = "NDVI") #calculate NDVI
  names(ndvi) <- i #rename layer
  if(!exists("ndvi_st")){
    ndvi_st <- ndvi
  } else {
    ndvi_st <- stack(ndvi_st, ndvi) #stack layers
  }
}
rm(ndvi);rm(tmp_stack);rm(i);rm(files_planetscope)

## crop raster stack to single trees and create long format dataframe 
dates <- unique(substr(names(ndvi_st),2,9)) #dates of planetscope images

ndvi_long <- data.frame(tree_id = NULL, date = NULL, values = NULL)
for(tree in as.character(unique(trees$id))){
  las1 <- readLAS(paste0("C:/Users/hhans/Documents/lidar_crowns/data/trees/",tree,".las"))
  las_shp <- lidR::as.spatial(las1)
  ndvi_st_crop <- crop(ndvi_st, las_shp)
  for(date in dates){
    for(i in which(substr(names(ndvi_st_crop), 2, 9) == date)){
      values <- ndvi_st_crop[[i]]@data@values
      ndvi_long <- rbind(ndvi_long, data.frame(tree_id=rep(tree, length(values)),
                                                 date=rep(as.Date(date, format="%Y%m%d"), length(values)),
                                                 values = values))
    }
  }
  
}
rm(ndvi_st_crop);rm(i);rm(date);rm(dates)
ndvi_long <- merge(ndvi_long, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F) #merge with phenoclasses

#################################################################
## create boxplots; timeseries of NDVI; grouped by phenophase ##
###############################################################
## boxplots - single tree
for(tree in unique(ndvi_long$tree_id)){
  out <- ndvi_long %>% filter(tree_id %in% tree) %>% 
    ggplot(aes(x=date, y=values, group=date, fill=as.factor(budburst_perc))) +
    geom_boxplot() +
    scale_fill_manual(values= cbPalette) +
    scale_x_date(date_breaks = "1 week") +
    xlab("Date") +
    ylab("NDVI") +
    labs(fill="buds bursted (in %)") +
    theme_light()
  ggsave(paste0("out/planetscape/single_tree_ndvi_timelines/",tree, "_ndvi_boxplot.png"), out, height = 10, width = 15)
}

  
## boxplots - all trees
out <- ndvi_long %>% 
  ggplot(aes(x=date, y=values, group=tree_id, fill=as.factor(budburst_perc))) +
  geom_boxplot(aes(group=date)) +
  facet_grid(tree_id~., scales = "free_y") +
  scale_fill_manual(values= cbPalette) +
  scale_x_date(date_breaks = "1 week") +
  xlab("Date") +
  ylab("NDVI") +
  labs(fill="buds bursted (in %)") +
  theme_light()
ggsave(paste0("out/planetscape/all_trees_ndvi_timeline/all_trees_ndvi_boxplot.png"), out, height = 10, width = 15)


#################################################################
## create plots; timeseries of NDVI; grouped by phenophase ##
###############################################################
## summarize data to tree, phenophase and date and calculate standard deviation
# function found here: http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

ndvi_sum <- data_summary(ndvi_long, varname="values", 
                    groupnames=c("tree_id", "date", "budburst", "budburst_perc"))

## ggplots points errorbars - single trees
for(tree in unique(ndvi_long$tree_id)){
  out <- ndvi_sum %>% filter(tree_id %in% tree) %>% 
    ggplot(aes(x=date, y=values, group=date, color=as.factor(budburst_perc))) +
    geom_point(size = 3)+
    geom_line() +
    geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=1, size = 1, position=position_dodge(0.05)) +
    scale_color_manual(values= cbPalette) +
    scale_x_date(date_breaks = "1 week") +
    xlab("Date") +
    ylab("NDVI") +
    labs(color="buds bursted (in %)") +
    theme_light()
  ggsave(paste0("out/planetscape/single_tree_ndvi_timelines/",tree, "_ndvi_points_error.png"), out, height = 10, width = 15)
}

## ggplots points errorbars - all trees
out <- ndvi_sum %>% 
  ggplot(aes(x=date, y=values, group=tree_id, color=as.factor(budburst_perc))) +
  geom_point(size = 3)+
  geom_line() +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=1, size = 1, position=position_dodge(0.05)) +
  facet_grid(tree_id~., scales = "free_y") +
  scale_color_manual(values= cbPalette) +
  scale_x_date(date_breaks = "1 week") +
  xlab("Date") +
  ylab("NDVI") +
  labs(color="buds bursted (in %)") +
  theme_light()
ggsave(paste0("out/planetscape/all_trees_ndvi_timeline/all_trees_ndvi_points_error.png"), out, height = 10, width = 15)
