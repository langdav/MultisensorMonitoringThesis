rm(list = ls())
library(plyr);library(dplyr);library(ggplot2);library(viridis);library(tidyr)

## create colorblind friendle palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## load data
# phenoclasses
phenoclasses <- read.csv("data/budburst_data/budburst_long.csv")
phenoclasses$date <- as.Date(phenoclasses$date)

# trees
load("data/trees_all.RData")
trees <- st_transform(trees, 25832)
trees$tree_id <- as.character(trees$tree_id)
#trees <- trees %>% filter(tree_id %in% c("mof_cst_00001","mof_cst_00003","mof_cst_00006","mof_cst_00013","mof_cst_00032","mof_cst_00036","mof_cst_00050", "BSF_1"))
trees <- trees[c(1:50),] #reduce to trees with a budburst record

# indices; Sentinel-1
sen1_indices <- read.csv("out/sentinel1/indicesi_long_format_phenoclasses_sentinel1.csv")
sen1_indices$date <- as.Date(sen1_indices$date)

# NDVI; Sentinel-2
sen2_ndvi <- read.csv("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.csv")
sen2_ndvi$date <- as.Date(sen2_ndvi$date)

# NDVI; planetscope
plan_ndvi <- read.csv("out/planetscope/ndvi_long_format_phenoclasses_planetscape.csv")
plan_ndvi$date <- as.Date(plan_ndvi$date)
names(plan_ndvi) <- c("tree_id","date","ndvi","budburst","budburst_perc")

# NDVI; orthomosaic,
#ortho_ndvi <- read.csv("out/orthomosaic/ndvi_per_tree_20210504.rds")

###### compare single trees in various products
#test <- sen2_ndvi %>% filter(tree_id %in% unique(trees$tree_id)[1]) %>% dplyr::select(ndvi)
sen2_tmp <- sen2_ndvi %>% filter(tree_id %in% unique(trees$tree_id)[1])
sen2_tmp$platform <- rep("sentinel2", nrow(sen2_tmp))
plan_tmp <- rbind(tmp, plan_ndvi %>% filter(tree_id %in% unique(trees$tree_id)[1]))
plan_tmp$platform <- rep("planetscope", nrow(plan_tmp))

platforms_tmp <- rbind(sen2_tmp, plan_tmp)


#########################
platforms_tmp %>%
  ggplot(aes(x=date, y=ndvi, color=platform)) +
  geom_boxplot(size = 1) +
  theme_light() +
  xlab("Date") +
  labs(fill = "Budburst Phase") +
  labs(color = "Tree ID")


platforms_tmp2 <- platforms_tmp %>% filter(date %in% as.Date(c("2021-03-20","2021-03-30","2021-05-09")))

platforms_tmp2 %>% 
  ggplot(aes(x=date, y=ndvi, color=platform)) +
  geom_point()

mean(platforms_tmp$)