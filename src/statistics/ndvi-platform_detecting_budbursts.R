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




