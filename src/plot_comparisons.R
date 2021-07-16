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



