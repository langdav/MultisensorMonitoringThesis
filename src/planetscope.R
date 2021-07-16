library(plyr);library(dplyr);library(ggplot2);library(viridis);library(tidyr)
library(raster);library(rgeos);library(rgdal);library(lidR);library(sf)

## load data
# phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

# trees
# trees <- readRDS("data/trees.RDS")
load("data/trees_all.RData")
trees <- st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
#trees <- trees %>% filter(tree_id %in% c("mof_cst_00001","mof_cst_00003","mof_cst_00006","mof_cst_00013","mof_cst_00032","mof_cst_00036","mof_cst_00050", "BSF_1"))
trees <- trees[c(1:50),] #reduce to trees with a budburst record

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
for(tree in as.character(unique(trees$tree_id))){
  
  las_shp <- rgdal::readOGR(paste0("data/single_tree_shapefiles/",tree,".gpkg"))
  crs(las_shp) <- CRS("+proj=longlat +datum=WGS84")
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
ndvi_long <- merge(ndvi_long, phenoclasses, by = c("tree_id","date"), all.x = F, all.y = F) #merge with phenoclasses
write.csv(ndvi_long, "out/planetscape/ndvi_long_format_phenoclasses_planetscape", row.names = F)

#################################################################
## create boxplots; timeseries of NDVI; grouped by phenophase ##
###############################################################
## boxplots - single tree
ndvi_long <- read.csv("data/ndvi_long.csv")
ndvi_long$date <- as.Date(ndvi_long$date)

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
  ggsave(paste0("out/planetscape/plots/single_tree_ndvi_timelines/",tree, "_ndvi_boxplot.png"), out, height = 10, width = 15)
}

  
## boxplots - all trees
ndvi_long <- read.csv("data/ndvi_long.csv")
ndvi_long$date <- as.Date(ndvi_long$date)

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
ggsave(paste0("out/planetscape/plots/all_trees_ndvi_timeline/all_trees_ndvi_boxplot.png"), out, height = 10, width = 15)


#################################################################
## create plots; timeseries of NDVI; grouped by phenophase ##
###############################################################
## summarize data to tree, phenophase and date and calculate standard deviation
# function found here: http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}

ndvi_long <- read.csv("data/ndvi_long.csv")
ndvi_long$date <- as.Date(ndvi_long$date)

ndvi_sum <- data_summary(ndvi_long, varname="values", 
                    groupnames=c("tree_id", "date", "budburst", "budburst_perc"))


## ggplots points errorbars - single trees
for(tree in unique(ndvi_long$tree_id)){
  out <- ndvi_sum %>% filter(tree_id %in% tree) %>% 
    ggplot(aes(x=date, y=values)) +
    geom_point(size = 3, aes(color=as.factor(budburst_perc)))+
    geom_errorbar(aes(ymin=values-sd, ymax=values+sd, color=as.factor(budburst_perc)), width=1, size = 1, position=position_dodge(0.05)) +
    geom_smooth(data = ndvi_long[which(ndvi_long$tree_id == tree),], aes(y=values), col = "black", size = .8) +
    scale_color_manual(values= cbPalette) +
    scale_x_date(date_breaks = "1 week") +
    xlab("Date") +
    ylab("NDVI") +
    labs(color="buds bursted (in %)") +
    theme_light()
  ggsave(paste0("out/planetscape/plots/single_tree_ndvi_timelines/",tree, "_ndvi_points_error.png"), out, height = 10, width = 15)
}


## ggplots points errorbars - all trees
out <- ndvi_sum %>% 
  ggplot(aes(x=date, y=values, group=tree_id, color=as.factor(budburst_perc))) +
  geom_point(size = 3)+
  geom_line() +
  geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=1, size = 1, position=position_dodge(0.05)) +
  #facet_grid(tree_id~., scales = "free_y") +
  #scale_color_manual(values= cbPalette) +
  scale_color_viridis(discrete = T) +
  scale_x_date(date_breaks = "1 week") +
  xlab("Date") +
  ylab("NDVI") +
  labs(color="buds bursted (in %)") +
  theme_light()
ggsave(paste0("out/planetscape/plots/all_trees_ndvi_timeline/all_trees_ndvi_points_error.png"), out, height = 10, width = 15)

######################################################
## NDVI points + fitted (loess, formula y + x ); per tree
ndvi_long %>% 
  ggplot(aes(date, values)) +
  geom_smooth(method = "loess", formula = y~x) +
  geom_point(aes(color=as.factor(budburst))) +
  #scale_color_manual(name = "Budburst", values= cbPalette) +
  scale_color_viridis(name = "Budburst", discrete = T) +
  scale_x_date(date_breaks = "1 week") +
  facet_grid(tree_id~., scales = "free_y") +
  xlab("Date") +
  ylab("NDVI") +
  labs(color="buds bursted (in %)") +
  theme_light()
ggsave(paste0("out/planetscape/plots/fitted_ndvi_budburst_phases_per_tree.png"), last_plot(), height = 10, width = 15)

## NDVI points + fitted (loess, formula y + x ); all trees; 26400 data points!
ndvi_long %>% 
  ggplot(aes(date, values)) +
  geom_smooth(method = "loess", formula = y~x, col = "black", size = .8) +
  geom_point(aes(color=as.factor(budburst))) +
  #scale_color_manual(name = "Budburst", values= cbPalette) +
  scale_color_viridis(name = "Budburst", discrete = T) +
  scale_x_date(date_breaks = "1 week") +
  xlab("Date") +
  ylab("NDVI") +
  labs(color="buds bursted (in %)") +
  theme_light()
ggsave(paste0("out/planetscape/plots/fitted_ndvi_budburst_phases_all_trees.png"), last_plot(), height = 10, width = 15)

## budburst percent per tree + fitted NDVI (loess, formula y + x)
ndvi_long %>% 
  ggplot(aes(date, budburst_perc)) +
  geom_line(aes(group = tree_id, color = tree_id), size = 1) +
  geom_smooth(aes(y = 100/max(values, na.rm = T)*values), method = "loess", col = "black", size = .8) +
  scale_x_date(date_breaks = "1 week") +
  scale_y_continuous(
    name = "Percentage of bursted buds",
    sec.axis = sec_axis(trans = ~.*(0.9237428/100),
                        name = "NDVI")
  ) +
  xlab("Date") +
  ylab("NDVI") +
  labs(color="Tree ID") +
  theme_light()
ggsave(paste0("out/planetscape/plots/fitted_ndvi_budburst_phases_per_tree.png"), last_plot(), height = 10, width = 15)

