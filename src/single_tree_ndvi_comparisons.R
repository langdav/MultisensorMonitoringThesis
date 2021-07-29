rm(list = ls())
library(plyr);library(dplyr);library(ggplot2);library(viridis)

# compare daily per tree NDVI

################################################################################
## Planetscope
load("out/planetscope/ndvi_long_format_phenoclasses_planetscope.RData")
planetscope <- ndvi_long_pheno_planetscope; rm(ndvi_long_pheno_planetscope)

# remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(planetscope$tree_id)){
  days <- unique(planetscope$date[which(planetscope$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(planetscope$ndvi[which(planetscope$tree_id == tree & 
                                                   planetscope$date == days[day])])$out) != 0){
      planetscope <- planetscope[-which(planetscope$tree_id == tree & 
                                          planetscope$date == days[day] & 
                                          planetscope$ndvi %in% c(boxplot.stats(planetscope$ndvi[which(planetscope$tree_id == tree & 
                                                                                                         planetscope$date == days[day])])$out)),]
    }
  }
}

# after removing outliers, calculate daily mean NDVI values
planetscope_daily_mean <- planetscope %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

################################################################################
## Sentinel-2
load("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.RData")
sentinel2 <- ndvi_long_pheno_sentinel2; rm(ndvi_long_pheno_sentinel2)
sentinel2$date <- as.Date(sentinel2$date)
sentinel2 <- rename(sentinel2, ndvi_mean = ndvi)


################################################################################
## TreeTalker
load("out/treetalker/ndvi_long_format_phenoclasses_treetalker.RData")
treetalker <- ndvi_long_format_phenoclasses_treetalker; rm(ndvi_long_format_phenoclasses_treetalker)

# remove outliers within each trees single-day-NDVI-values, based on the IQR
for(tree in unique(treetalker$tree_id)){
  days <- unique(treetalker$date[which(treetalker$tree_id == tree)])
  for(day in 1:length(days)){
    if(length(boxplot.stats(treetalker$ndvi[which(treetalker$tree_id == tree & 
                                                   treetalker$date == days[day])])$out) != 0){
      treetalker <- treetalker[-which(treetalker$tree_id == tree & 
                                          treetalker$date == days[day] & 
                                          treetalker$ndvi %in% c(boxplot.stats(treetalker$ndvi[which(treetalker$tree_id == tree & 
                                                                                                         treetalker$date == days[day])])$out)),]
    }
  }
}

# after removing outliers, calculate daily mean NDVI values
treetalker_daily_mean <- treetalker %>% 
  group_by(tree_id, date) %>% 
  summarize(ndvi_mean = mean(ndvi, na.rm = T),
            ndvi_sd = sd(ndvi, na.rm = T),
            budburst = unique(budburst),
            budburst_perc = unique(budburst_perc))

################################################################################
## Orthomosaic
load("out/orthomosaic/ndvi_mean_per_tree_long_format_phenoclasses_orthomosaic.RData")
orthomosaic <- ndvi_mean_pheno_orthomosaic; rm(ndvi_mean_pheno_orthomosaic)
orthomosaic <- rename(orthomosaic, ndvi_mean = ndvi)

################################################################################
## Plots

# all in one
planetscope_daily_mean$platform <- rep("planetscope", nrow(planetscope_daily_mean))
sentinel2$platform <- rep("sentinel2", nrow(sentinel2))
treetalker_daily_mean$platform <- rep("treetalker", nrow(treetalker_daily_mean))
orthomosaic$platform <- rep("orthomosaic", nrow(orthomosaic))

aio <- rbind(planetscope_daily_mean, sentinel2, treetalker_daily_mean)

aio %>% filter(tree_id %in% "mof_cst_00017") %>% 
  ggplot(aes(x=date, y=ndvi_mean, color=platform)) +
  geom_line()

sentinel2 %>% 
  ggplot(aes(x=date, y=ndvi, group=tree_id, color=as.factor(budburst_perc))) +
  geom_point(size = 3)+
  geom_line() +
  #facet_grid(tree_id~., scales = "free_y") +
  #scale_color_manual(values= cbPalette) +
  scale_color_viridis(discrete = T) +
  scale_x_date(date_breaks = "1 week") +
  xlab("Date") +
  ylab("NDVI") +
  labs(color="buds bursted (in %)") +
  theme_light()

################################################################################
## test for normal distribution 
data <- orthomosaic$ndvi_mean
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F, main = paste(colnames(data)))
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

shapiro.test(data)
#### non of the data is normal distributed -> non parametric tests -> Anova

################################################################################
## testing for group differences -> different platforms; different trees; mean NDVI values over whole period 
# for Anova: testing if residuals are normally distributed
tree <- unique(aio$tree_id)[5]
model <- lm(ndvi_mean~platform, data = aio[which(aio$tree_id == tree),])
car::qqPlot(rstandard(model)) # somewhat normally distributed residuals

# Anova for single tree
anova_plan <- aov(ndvi_mean~platform, data = aio[which(aio$tree_id == tree),])
summary(anova_plan) # p > 0.05 -> means DO NOT differ

# Anova for all trees
out <- NULL
for(tree in unique(aio$tree_id)){
  anova <- aov(ndvi_mean~platform, data = aio[which(aio$tree_id == tree),])
  pval <- summary(anova)[[1]][["Pr(>F)"]][1]
  ifelse(pval < 0.05, different <- T, different <- F)
  out <- rbind(out, data.frame(tree_id = tree,
                               p_value = pval,
                               means_differ = different))
}
out

################################################################################
## testing for group differences -> different platforms; different days; mean NDVI values of whole forest (all trees)
# for Anova: testing if residuals are normally distributed
date <- unique(aio$date[which(aio$platform == "sentinel2")])[3] # first date is only available for sentinel2 data
model <- lm(ndvi_mean~platform, data = aio[which(aio$date == date),])
car::qqPlot(rstandard(model)) # somewhat normally distributed residuals

# Anova for single tree
anova_plan <- aov(ndvi_mean~platform, data = aio[which(aio$date == date),])
summary(anova_plan) # p > 0.05 -> means DO NOT differ

# Anova for all trees
out <- NULL
for(date in unique(aio$date[which(aio$platform == "sentinel2")])[2:3]){ #sentinel has the least amounts of days with available data; so just compare those days
  anova <- aov(ndvi_mean~platform, data = aio[which(aio$date == date),])
  pval <- summary(anova)[[1]][["Pr(>F)"]][1]
  ifelse(pval < 0.05, different <- T, different <- F)
  out <- rbind(out, data.frame(tree_id = tree,
                               p_value = pval,
                               means_differ = different))
}
out

