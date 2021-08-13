################################################################################
## Test:
# can the budburst be detected within different platforms per-tree mean NDVI values?
## Results: 
# 
################################################################################

rm(list = ls())
library(plyr);library(dplyr);library(ggplot2);library(viridis);library(rstatix);library(ggpubr);library(car)

# load data ---------------------------------------
# -------------------------------------------------
# NDVI means per tree, day and platform
load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)

# NDVI means over all trees, per day and platform
load("out/all_in_one/aio_all_tree_means_per_platform_and_day.RData")

# =============================================================


# Planetscope ---------------------------------------
# ---------------------------------------------------
planetscope <- aio %>% filter(platform == "planetscope")
planetscope_all_tree_means <- aio_all_tree_means %>% filter(platform == "planetscope")

# per tree NDVI
planetscope %>% 
  ggplot(aes(x=date, y=ndvi_mean, color = as.factor(budburst))) +
  geom_smooth(method = loess)

# all trees mean NDVI
planetscope_all_tree_means %>% 
  ggplot(aes(x=date, y=ndvi_mean, color = as.factor(budburst))) +
  geom_smooth(method = loess)


# =============================================================

# Orthomosaic ---------------------------------------
# ---------------------------------------------------
orthomosaic <- aio %>% filter(platform == "orthomosaic")
orthomosaic_all_tree_means <- aio_all_tree_means %>% filter(platform == "orthomosaic")

# per tree NDVI
orthomosaic %>% 
  ggplot(aes(x=date, y=ndvi_mean, color = as.factor(budburst))) +
  geom_smooth(method = loess)

# all trees mean NDVI
orthomosaic_all_tree_means %>% 
  ggplot(aes(x=date, y=ndvi_mean, color = as.factor(budburst))) +
  geom_smooth(method = loess)


# =============================================================

# Sentinel-2 ---------------------------------------
# ---------------------------------------------------
sentinel2 <- aio %>% filter(platform == "sentinel2")
sentinel2_all_tree_means <- aio_all_tree_means %>% filter(platform == "sentinel2")

# fitting a mopdel does not make sense here, as there is only one image containing budburst = T
sentinel2 %>% 
  ggplot(aes(x=date, y=ndvi_mean, color = as.factor(budburst))) +
  geom_boxplot(aes(group = as.factor(budburst)))

# all trees mean NDVI
sentinel2_all_tree_means %>% 
  ggplot(aes(x=date, y=ndvi_mean, color = as.factor(budburst))) +
  geom_smooth(method = loess)

# =============================================================

############################################################################################
## per sentinel1-indice plots; comparison of all 5 trees; budburst phases as background ##
##########################################################################################
## add another coloumn needed for geom_rects
# aio$budtime <- NA
# for(platform in unique(aio$platform)){
#   prows <- which(aio$platform == platform)
#   for(i in 1:length(prows)){
#     aio$budtime[prows[i]] <- ifelse(i < length(prows),as.character(aio$date[c(prows[i+1])]),as.character(aio$date[c(prows[i])]))
#   }
# }
# aio$budtime <- as.Date(aio$budtime)

## create data frame for "budburst" and "no budburst" geom_rects
budburst <- read.csv("data/budburst_data/budburst_long.csv")
buddies_budburst_classes <- NULL
for(tree in unique(budburst$tree_id)){
  tmp <- budburst %>% filter(tree_id == tree)
  buddies_budburst_classes <- rbind(buddies_budburst_classes, data.frame(tree_id = tree,
                                       first_date = as.Date(head(tmp$date,1)),
                                       last_date = as.Date(tmp$date[min(which(tmp$budburst==T))]),
                                       budburst = F))
  buddies_budburst_classes <- rbind(buddies_budburst_classes, data.frame(tree_id = tree,
                                       first_date = as.Date(tmp$date[min(which(tmp$budburst==T))]),
                                       last_date = as.Date(tail(tmp$date,1)),
                                       budburst = T))
  
}


## create data frame for budburst percent geom_rects
buddies_budburst_percent <- NULL
for(tree in unique(budburst$tree_id)){
  tmp <- budburst %>% filter(tree_id == tree)
  percs <- unique(tmp$budburst_perc)
  for(perc in percs){
    if(perc == 0){
      buddies_budburst_percent <- rbind(buddies_budburst_percent, data.frame(tree_id = tree,
                                                                             first_date = as.Date(head(tmp$date,1)),
                                                                             last_date = as.Date(tmp$date[max(which(tmp$budburst_perc==0))]),
                                                                             budburst_perc = 0))
    }
    if(perc == 100){
      buddies_budburst_percent <- rbind(buddies_budburst_percent, data.frame(tree_id = tree,
                                                                             first_date = as.Date(tmp$date[min(which(tmp$budburst_perc==100))]),
                                                                             last_date = as.Date(tail(tmp$date,1)),
                                                                             budburst_perc = 100))
    }
    if(perc > 0 & perc < 100){
      buddies_budburst_percent <- rbind(buddies_budburst_percent, data.frame(tree_id = tree,
                                                                             first_date = as.Date(tmp$date[min(which(tmp$budburst_perc==perc))]),
                                                                             last_date = as.Date(tmp$date[max(which(tmp$budburst_perc==perc))]),
                                                                             budburst_perc = perc))
    }
  }
}

## plot
tree <- unique(aio$tree_id)[4]

ggplot() +
  geom_line(data = aio[which(aio$tree_id==tree),], aes(x=date, y=ndvi_mean, group=platform, color=platform), size = 1) +
  facet_grid(platform~., scales = "free_y") +
  scale_color_viridis(discrete = T) + 
  theme_light() +
  xlab("Date") +
  scale_x_date(date_breaks = "2 weeks") +
  ylab(names(aio)[i]) +
  labs(fill = "Budburst Phase") +
  labs(color = "Tree ID")+
  geom_rect(data = buddies_budburst_classes[which(buddies_budburst_classes$tree_id == tree),], aes(xmin=first_date, xmax=last_date, ymin = 0, ymax=1, fill = budburst), color = NA, alpha=0.5) +
  scale_fill_brewer(palette="YlGn")  #Greys = RColorBrewer Code; alternative (good) colorspaces: YlGn, Greens


## plot
tree <- unique(aio$tree_id)[34]

ggplot() +
  geom_line(data = aio[which(aio$tree_id==tree),], aes(x=date, y=ndvi_mean, group=platform), size = 1) +
  facet_grid(platform~., scales = "free_y") +
  scale_color_viridis(discrete = T) + 
  theme_light() +
  xlab("Date") +
  scale_x_date(date_breaks = "2 weeks") +
  ylab(names(aio)[i]) +
  labs(fill = "Budburst Phase") +
  labs(color = "Tree ID")+
  geom_rect(data = buddies_budburst_percent[which(buddies_budburst_percent$tree_id == tree),], aes(xmin=first_date, xmax=last_date, ymin = 0, ymax=1, fill = as.factor(budburst_perc)), color = NA, alpha=0.5) +
  scale_fill_brewer(palette="YlGn")  #Greys = RColorBrewer Code; alternative (good) colorspaces: YlGn, Greens


### what does this show?
# if we are only distinguishing between "budburst" and "no budburst", which is a prerequisite if we are NOT looking onto single tree level, 
# we are missing major changes in NDVI; 
# graphs show, that that large jumps in NDVI are taking place in distinct levels of budburst percentage
# to see those, we NEED to look onto trees on a single tree level

## add slopes for NDVI over time for each platform and each tree
slope <- NULL
highest_slope <- NULL
for(platform in unique(aio$platform)){
  for(tree in unique(aio$tree_id[which(aio$platform == platform)])){
    slope_tmp <- c(diff(aio$ndvi_mean[which(aio$platform == platform & aio$tree_id==tree)]) / diff(as.numeric(aio$date[which(aio$platform == platform & aio$tree_id==tree)])))
    highest_slope_tmp <- rep(F, length(slope_tmp))
    highest_slope_tmp[which(slope_tmp == max(slope_tmp))-1] <- T
    slope <- c(slope, c(NA,slope_tmp))
    highest_slope <- c(highest_slope, c(F,highest_slope_tmp))
  }
}

aio_slope <- cbind(aio, slope, highest_slope)

tree <- unique(aio_slope$tree_id)[2]

ggplot() +
  geom_line(data = aio_slope[which(aio_slope$tree_id==tree),], aes(x=date, y=ndvi_mean, group=platform, color = highest_slope), size = 1) +
  facet_grid(platform~., scales = "free_y") +
  scale_color_discrete() +
  theme_light() +
  xlab("Date") +
  scale_x_date(date_breaks = "2 weeks") +
  ylab(names(aio_slope)[i]) +
  labs(fill = "Budburst Phase") +
  labs(color = "Highest Slope")+
  geom_rect(data = buddies_budburst_percent[which(buddies_budburst_percent$tree_id == tree),], aes(xmin=first_date, xmax=last_date, ymin = 0, ymax=1, fill = as.factor(budburst_perc)), color = NA, alpha=0.5) +
  scale_fill_brewer(palette="Greys")  #Greys = RColorBrewer Code; alternative (good) colorspaces: YlGn, Greens

aio_slope[which(aio_slope$tree_id==tree & aio$platform=="planetscope"),]
