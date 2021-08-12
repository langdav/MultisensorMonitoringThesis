################################################################################
## Test:
# checking for group differences
# groups: platforms "orthomosaic", "planetscope" and "sentinel2"
# values: per-tree daily mean NDVI values
# other stuff: comparing NDVI values within classes "budburst" & "no budburst"
## Results: 
# 
################################################################################

rm(list = ls())
library(plyr);library(dplyr);library(ggplot2);library(viridis);library(rstatix);library(ggpubr);library(car)


# load all-in-one daily means data frame ---------------------------------------
# ------------------------------------------------------------------------------
load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
# =============================================================


# outlier detection and removal ---------------------------------------
# ---------------------------------------------------------------------
## class: "budburst"
# platform: Planetscope
boxplot(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == T)]) #no outliers

# platform: Orthomosaic
boxplot(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == T)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == T)])$out
aio <- aio[-which(aio$platform == "orthomosaic" &
                    aio$budburst == T & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == T)])$out),]

# platform: Sentinel-2
boxplot(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)])$out
aio <- aio[-which(aio$platform == "sentinel2" &
                    aio$budburst == T & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)])$out),]

## class "no budburst"
# platform: Planetscope
boxplot(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == F)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == F)])$out
aio <- aio[-which(aio$platform == "planetscope" &
                    aio$budburst == F & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == F)])$out),]

# platform: Orthomosaic
boxplot(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == F)]) # no outliers

# platform: Sentinel-2
boxplot(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == F)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == F)])$out
aio <- aio[-which(aio$platform == "sentinel2" &
                    aio$budburst == F & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == F)])$out),]

# =============================================================

# check for normal distribution ---------------------------------------
# ---------------------------------------------------------------------
#data <- aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == T)] #not normally distributed
#data <- aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == T)] #normally distributed
#data <- aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)] #normally distributed
#data <- aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == F)] #not normally distributed
#data <- aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == F)] #not normally distributed
#data <- aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)] #normally distributed
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F, main = paste(colnames(data)))
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

rstatix::shapiro_test(data) # not all data are normally distributed -> Anova

### Anova prerequisite: testing if residuals are normally distributed
## test for normal distribution of residuals within all groups (fast and easy)
#model  <- lm(ndvi_mean ~ platform, data = aio[which(aio$budburst == F),]) #p < 0.5; residuals not normally distributed
#model  <- lm(ndvi_mean ~ platform, data = aio[which(aio$budburst == T),]) #p < 0.5; residuals not normally distributed
ggpubr::ggqqplot(residuals(model))
rstatix::shapiro_test(residuals(model))
# as the sample sizes per group are big enough (1790 for budburst = T and 2215 for budburst = F), an Anova may be used although the residuals are not normally distributed

## test for normal distribution for each group individually
ggpubr::ggqqplot(aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),], "ndvi_mean", facet.by = "platform")

### Anova prerequisite: check for variance homogeneity
car::leveneTest(ndvi_mean~platform, data = aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),]) # p < 0.05; no Variance homogeneity -> WELCH-ANOVA

# =============================================================


### ANOVA; Welch-Anova, as variances are not equal ---------------------------------------
# ----------------------------------------------------------------------------------------
## group "budburst"
aio %>% 
  filter(budburst == T) %>% 
  filter(platform %in% c("orthomosaic", "planetscope", "sentinel2")) %>% 
  group_by() %>% 
  select(ndvi_mean, platform) %>% 
  rstatix::welch_anova_test(ndvi_mean ~ platform)
# p < 0.05 -> groups are significantly different

## Tukey post-hoc test
aio %>% 
  filter(budburst == T) %>% 
  filter(platform %in% c("orthomosaic", "planetscope", "sentinel2")) %>% 
  group_by() %>% 
  select(ndvi_mean, platform) %>% 
  rstatix::tukey_hsd(ndvi_mean ~ platform)
# orthomosaic - planetscope; not signif.
# orthomosaic - sentinel2; signif.
# planetscope - sentinel2; signif.

## group " no budburst"
aio %>% 
  filter(budburst == F) %>% 
  filter(platform %in% c("orthomosaic", "planetscope", "sentinel2")) %>% 
  group_by() %>% 
  select(ndvi_mean, platform) %>% 
  rstatix::welch_anova_test(ndvi_mean ~ platform)
# p < 0.05 -> groups are significantly different

## Tukey post-hoc test
aio %>% 
  filter(budburst == F) %>% 
  filter(platform %in% c("orthomosaic", "planetscope", "sentinel2")) %>% 
  group_by() %>% 
  select(ndvi_mean, platform) %>% 
  rstatix::tukey_hsd(ndvi_mean ~ platform)
# orthomosaic - planetscope; signif.
# orthomosaic - sentinel2; signif.
# planetscope - sentinel2; signif.

# =============================================================
