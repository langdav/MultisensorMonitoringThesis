#author: David Langenohl
#last modified: 17.08.2021
#description: calculate NDVI values from Planetscope images
#NOTE: blue = 1, green = 2, red = 3, nir = 4

rm(list = ls())
library(plyr);library(dplyr);library(ggplot2);library(viridis);library(tidyr)
library(raster);library(rgeos);library(rgdal);library(lidR);library(sf)

## list available files
files_planetscope <- list.files("data/satellite_data/planetscope/",pattern = "SR_clip.tif")

## calculate ndvi for all trees
ndvi_all <- data.frame(tree_id = NULL, date = NULL, ndvi = NULL)

for(file in files_planetscope){
  tmp_stack <- stack(paste0("data/satellite_data/planetscope/", file))
  
  load("data/trees.RData")
  trees$tree_id <- as.character(trees$tree_id)
  trees <- trees[c(1:50),] #reduce to trees with a budburst record
  
  for(i in 1:nrow(trees)){
    cat("Processing", trees$tree_id[i], "in file", file,"\n")
    
    #load single tree shapefile
    single_tree_sf <- sf::read_sf(paste0("data/single_tree_shapefiles/",trees$tree_id[i],".gpkg"))
    
    # convert Spatial Points to Polygon Shapefile
    single_tree_sf <- single_tree_sf$geom %>% 
      st_union() %>% 
      st_convex_hull()
    single_tree_sf <- st_as_sf(single_tree_sf)
    
    #extract green, red and NIR values for single tree
    green <- raster::extract(tmp_stack[[2]], single_tree_sf, na.rm = T)[[1]]
    red <- raster::extract(tmp_stack[[3]], single_tree_sf, na.rm = T)[[1]]
    nir <- raster::extract(tmp_stack[[4]], single_tree_sf, na.rm = T)[[1]]
    
    #calculate NDVI = (nir-red)/(nir+red)
    ndvi <- (nir-red)/(nir+red)
    
    #create output data frame
    ndvi_all <- rbind(ndvi_all, data.frame(tree_id=rep(trees$tree_id[i], length(ndvi)),
                                           date=rep(as.Date(substr(file, 1, 8), "%Y%m%d"), length(ndvi)),
                                           ndvi = ndvi))
    
  }
}

#merge with phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

planetscope_ndvi_all <- merge(ndvi_all, phenoclasses, by = c("tree_id","date"), all.x = T, all.y = F)

#save resulting data frame
save(planetscope_ndvi_all, file = "out/planetscope/ndvi_all_with_phenoclasses_planetscope.RData")

# remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(planetscope_ndvi_all$tree_id)){
  cat("Processing", tree, "\n")
  days <- unique(planetscope_ndvi_all$date[which(planetscope_ndvi_all$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(planetscope_ndvi_all$ndvi[which(planetscope_ndvi_all$tree_id == tree & 
                                                   planetscope_ndvi_all$date == days[day])])$out) != 0){
      planetscope_ndvi_all <- planetscope_ndvi_all[-which(planetscope_ndvi_all$tree_id == tree & 
                                          planetscope_ndvi_all$date == days[day] & 
                                          planetscope_ndvi_all$ndvi %in% c(boxplot.stats(planetscope_ndvi_all$ndvi[which(planetscope_ndvi_all$tree_id == tree & 
                                                                                                         planetscope_ndvi_all$date == days[day])])$out)),]
    }
  }
}

#save resulting data frame
save(planetscope_ndvi_all, file = "out/planetscope/outlier_free_ndvi_all_with_phenoclasses_planetscope.RData")

# after removing outliers, calculate daily mean NDVI values
planetscope_ndvi_mean_per_tree <- planetscope_ndvi_all %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

# save resulting data frame
save(planetscope_ndvi_mean_per_tree, file = "out/planetscope/outlier_free_daily_ndvi_mean_per_tree_with_phenoclasses_planetscope.RData")



#####################################################################################################################################
#################################################################
## create boxplots; timeseries of NDVI; grouped by phenophase ##
###############################################################
## create colorblind friendle palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## boxplots - single tree
load("out/planetscope/ndvi_long_format_phenoclasses_planetscope.RData")
ndvi_long <- ndvi_long_pheno_planetscope;rm(ndvi_long_pheno_planetscope)

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
  ggsave(paste0("out/planetscope/plots/single_tree_ndvi_timelines/",tree, "_ndvi_boxplot.png"), out, height = 10, width = 15)
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
ggsave(paste0("out/planetscope/plots/all_trees_ndvi_timeline/all_trees_ndvi_boxplot.png"), out, height = 10, width = 15)


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

ndvi_long <- read.csv("out/planetscope/ndvi_long_format_phenoclasses_planetscape")
ndvi_long$date <- as.Date(ndvi_long$date)

ndvi_sum <- data_summary(ndvi_long, varname="ndvi", 
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
  ggsave(paste0("out/planetscope/plots/single_tree_ndvi_timelines/",tree, "_ndvi_points_error.png"), out, height = 10, width = 15)
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
ggsave(paste0("out/planetscope/plots/all_trees_ndvi_timeline/all_trees_ndvi_points_error.png"), out, height = 10, width = 15)

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
ggsave(paste0("out/planetscope/plots/fitted_ndvi_budburst_phases_per_tree.png"), last_plot(), height = 10, width = 15)

## NDVI points + fitted (loess, formula y + x ); all trees; 23000 data points!
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
ggsave(paste0("out/planetscope/plots/fitted_ndvi_budburst_phases_all_trees.png"), last_plot(), height = 10, width = 15)

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
ggsave(paste0("out/planetscope/plots/fitted_ndvi_budburst_phases_per_tree.png"), last_plot(), height = 10, width = 15)

