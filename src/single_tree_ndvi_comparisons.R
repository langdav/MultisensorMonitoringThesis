################################################################################
# Results:
# 1) all per tree daily NDVI (mean) values as well as all per platform all trees (whole forest) NDVI values are not normally distributed -> Anova
# 2) group differences -> different platforms; different trees; mean NDVI values over whole period
#### for most trees, the whole-period-NDVI-means differ signficantly
#### MAYBE DO A MORE IN-DEPTH TEST AND COMPARE PLATFORMS ONE BY ONE
# 3) group differences -> different platforms; different days; mean NDVI values of whole forest (all trees)
#### could not compare all platforms at once, as there is no date that is available in all platforms -> pairwise comparison
#### only compared remote sensing platforms: Sentinel-2, Planetscope, Orthomosaics
#### Sentinel-2 - Planetscope: NDVI-means of Sentinel2 and Planetscope differ on 2021-03-30; do not differ on 2021-05-09
#### Sentinel-2 - Orthomosaics: no common dates
#### Planetscope - Orthomosaics: NDVI-means of Orthomosaics and Planetscope differ on all three dates
################################################################################

rm(list = ls())
library(plyr);library(dplyr);library(ggplot2);library(viridis)

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
# per-tree mean NDVI values
load("out/orthomosaic/ndvi_mean_per_tree_long_format_phenoclasses_orthomosaic.RData")
orthomosaic_daily_means <- ndvi_mean_pheno_orthomosaic; rm(ndvi_mean_pheno_orthomosaic)
orthomosaic_daily_means <- rename(orthomosaic_daily_means, ndvi_mean = ndvi)

# all NDVI values (multiple values per tree)
load("out/orthomosaic/outlier_free_ndvi_mean_per_tree_long_format_phenoclasses_orthomosaic.RData")
# load("out/orthomosaic/ndvi_long_format_phenoclasses_orthomosaic.RData") #only load if necessary, as resulting dataframe is huge and contains 8553222 rows
# orthomosaic <- ndvi_long_pheno_orthomosaic; rm(ndvi_long_pheno_orthomosaic)
# 
# # remove outliers within each trees single-day-NDVI-values, based on the IQR
# for(tree in unique(orthomosaic$tree_id)){
#   days <- unique(orthomosaic$date[which(orthomosaic$tree_id == tree)])
#   for(day in 1:length(days)){
#     if(length(boxplot.stats(orthomosaic$ndvi[which(orthomosaic$tree_id == tree & 
#                                                    orthomosaic$date == days[day])])$out) != 0){
#       orthomosaic <- orthomosaic[-which(orthomosaic$tree_id == tree & 
#                                         orthomosaic$date == days[day] & 
#                                         orthomosaic$ndvi %in% c(boxplot.stats(orthomosaic$ndvi[which(orthomosaic$tree_id == tree & 
#                                                                                                        orthomosaic$date == days[day])])$out)),]
#     }
#   }
# }
# 
# # as outlier removal takes rather long, due to high number of entries, save outlier-free data frame and directly load that one the next times
# save(orthomosaic, file = "out/orthomosaic/outlier_free_ndvi_mean_per_tree_long_format_phenoclasses_orthomosaic.RData")

 ################################################################################
## Create all-in-one dataframe
planetscope_daily_mean$platform <- rep("planetscope", nrow(planetscope_daily_mean))
sentinel2$platform <- rep("sentinel2", nrow(sentinel2))
treetalker_daily_mean$platform <- rep("treetalker", nrow(treetalker_daily_mean))
orthomosaic_daily_means$platform <- rep("orthomosaic", nrow(orthomosaic_daily_means))

aio <- rbind(planetscope_daily_mean, sentinel2, treetalker_daily_mean,orthomosaic_daily_means)

################################################################################
## Plots

aio %>% filter(tree_id %in% "mof_cst_00017") %>% 
  ggplot(aes(x=date, y=ndvi_mean, color=platform)) +
  geom_line()

sentinel2 %>% 
  ggplot(aes(x=date, y=ndvi_mean, group=tree_id, color=as.factor(budburst_perc))) +
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
################################################################################
################################################################################
# 2) group differences -> different platforms; different trees; mean NDVI values over whole period

# test for normal distribution 
data <- orthomosaic_daily_means$ndvi_mean
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F, main = paste(colnames(data)))
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

shapiro.test(data)
#### non of the data is normal distributed -> non parametric tests -> Anova

# for Anova: testing if residuals are normally distributed
tree <- unique(aio$tree_id)[1]
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
################################################################################
################################################################################
# 3) group differences -> different platforms; different days; mean NDVI values of whole forest (all trees)

# seeing, which dates are available for which platforms
unique(orthomosaic_daily_means$date)
unique(planetscope_daily_mean$date)
unique(sentinel2$date)
unique(treetalker_daily_mean$date)

# comparison table
all_dates <- unique(aio$date)
com_dates <- data.frame(date = all_dates,
                        orthomosaic = F,
                        planetscope = F,
                        sentinel2 = F,
                        treetalker = F,
                        available_in_all_platforms = F,
                        number_of_available_platforms = NA)

com_dates$orthomosaic <- replace(com_dates$orthomosaic, which(com_dates$date %in% unique(orthomosaic_daily_means$date)), T)
com_dates$planetscope <- replace(com_dates$planetscope, which(com_dates$date %in% unique(planetscope_daily_mean$date)), T)
com_dates$sentinel2 <- replace(com_dates$sentinel2, which(com_dates$date %in% unique(sentinel2$date)), T)
com_dates$treetalker <- replace(com_dates$treetalker, which(com_dates$date %in% unique(treetalker_daily_mean$date)), T)

for(i in 1:nrow(com_dates)){
  if(all(com_dates[i,] == T)){
    com_dates$available_in_all_platforms[i] <- T
  }
  com_dates$number_of_available_platforms[i] <- length(which(com_dates[i,c(2:5)] == T))
}

# any dates, which are available for all platforms?
any(com_dates$available_in_all_platforms == T) #no

#############################################
## comparing all platforms, does not make sense here, as there is no date, that is available in all platforms
## But: I can do pairwise comparisons
#############################################
################################################################################
## Orthomosaic - Planetscope ################

# find dates, that both platforms have data for
available_dates <- com_dates$date[which(com_dates$orthomosaic == T & com_dates$planetscope == T)]

# test for normal distribution 
#data <- aio$ndvi_mean[which(aio$date %in% available_dates & aio$platform == "planetscope")]
data <- aio$ndvi_mean[which(aio$date %in% available_dates & aio$platform == "orthomosaic")]
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F, main = paste(colnames(data)))
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

shapiro.test(data) # non of the data is normally distributed -> non parametric tests -> Anova

# for Anova: testing if residuals are normally distributed; manually for every date
date <- available_dates[3] # first date is only available for sentinel2 data
model <- lm(ndvi_mean~platform, data = aio[which(aio$date == date & aio$platform %in% c("orthomosaic", "planetscope")),])
car::qqPlot(rstandard(model))

## residuals are normally distributed; Anova can be used
# Anova for all trees
out <- NULL
for(date in available_dates){ #sentinel has the least amounts of days with available data; so just compare those days
  anova <- aov(ndvi_mean~platform, data = aio[which(aio$date == date & aio$platform %in% c("orthomosaic", "planetscope")),])
  pval <- summary(anova)[[1]][["Pr(>F)"]][1]
  ifelse(pval < 0.05, different <- T, different <- F)
  out <- rbind(out, data.frame(date = unique(aio$date[which(aio$date == date)]),
                               p_value = pval,
                               means_differ = different))
}
out #NDVI-means of Orthomosaics and Planetscope differ on all three dates


################################################################################
## Orthomosaic - Sentinel2 ################

# find dates, that both platforms have data for
available_dates <- com_dates$date[which(com_dates$orthomosaic == T & com_dates$sentinel2 == T)] #no common dates


################################################################################
## Sentinel2 - Planetscope ################

# find dates, that both platforms have data for
available_dates <- com_dates$date[which(com_dates$sentinel2 == T & com_dates$planetscope == T)]

# test for normal distribution 
#data <- aio$ndvi_mean[which(aio$date %in% available_dates & aio$platform == "planetscope")]
data <- aio$ndvi_mean[which(aio$date %in% available_dates & aio$platform == "sentinel2")]
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F, main = paste(colnames(data)))
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

shapiro.test(data) # non of the data is normally distributed -> non parametric tests -> Anova

# for Anova: testing if residuals are normally distributed; manually for every date
date <- available_dates[2] # first date is only available for sentinel2 data
model <- lm(ndvi_mean~platform, data = aio[which(aio$date == date & aio$platform %in% c("sentinel2", "planetscope")),])
car::qqPlot(rstandard(model))

## residuals are normally distributed; Anova can be used
# Anova for all trees
out <- NULL
for(date in available_dates){ #sentinel has the least amounts of days with available data; so just compare those days
  anova <- aov(ndvi_mean~platform, data = aio[which(aio$date == date & aio$platform %in% c("sentinel2", "planetscope")),])
  pval <- summary(anova)[[1]][["Pr(>F)"]][1]
  ifelse(pval < 0.05, different <- T, different <- F)
  out <- rbind(out, data.frame(date = unique(aio$date[which(aio$date == date)]),
                               p_value = pval,
                               means_differ = different))
}
out #NDVI-means of Sentinel2 and Planetscope differ on 2021-03-30; do not differ on 2021-05-09


################################################################################
################################################################################
################################################################################
# 4) group differences -> different platforms; NDVI values within budburst classes "budburst" & "no budburst"

## outlier detection and removal
# budburst
boxplot(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == T)]) #no outliers

boxplot(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == T)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == T)])$out
aio <- aio[-which(aio$platform == "orthomosaic" &
                    aio$budburst == T & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == T)])$out),]

boxplot(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)])$out
aio <- aio[-which(aio$platform == "sentinel2" &
                    aio$budburst == T & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == T)])$out),]

# no budburst
boxplot(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == F)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == F)])$out
aio <- aio[-which(aio$platform == "planetscope" &
                    aio$budburst == F & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "planetscope" & aio$budburst == F)])$out),]

boxplot(aio$ndvi_mean[which(aio$platform == "orthomosaic" & aio$budburst == F)]) # no outliers

boxplot(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == F)]) # outliers
boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == F)])$out
aio <- aio[-which(aio$platform == "sentinel2" &
                    aio$budburst == F & 
                    aio$ndvi_mean %in% boxplot.stats(aio$ndvi_mean[which(aio$platform == "sentinel2" & aio$budburst == F)])$out),]

## check for normal distribution
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

shapiro.test(data) # not all data are normally distributed -> Anova

## for Anova: testing if residuals are normally distributed
#model <- lm(ndvi_mean~date, data = aio[which(aio$platform == "planetscope" & aio$budburst == T),]) # normally distributed
#model <- lm(ndvi_mean~date, data = aio[which(aio$platform == "orthomosaic" & aio$budburst == T),]) # normally distributed
#model <- lm(ndvi_mean~date, data = aio[which(aio$platform == "sentinel2" & aio$budburst == T),]) # normally distributed
#model <- lm(ndvi_mean~date, data = aio[which(aio$platform == "planetscope" & aio$budburst == F),]) # normally distributed
#model <- lm(ndvi_mean~date, data = aio[which(aio$platform == "orthomosaic" & aio$budburst == F),]) # normally distributed
#model <- lm(ndvi_mean~date, data = aio[which(aio$platform == "sentinel2" & aio$budburst == F),]) # normally distributed
model <- lm(ndvi_mean~platform, data = aio[which(aio$budburst == F),]) # normally distributed
car::qqPlot(model$residuals)

## check for variance homogeneity (Anova prerequisite)
library(car)
leveneTest(ndvi_mean~platform, data = aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),]) 
# p < 0.05; no Variance homogeneity -> WELCH-ANOVA

## check for normally distributed residuals within groups (Anova prerequisite)
# test for normal distribution of residuals within all groups (fast and easy)
model  <- lm(ndvi_mean ~ platform, data = aio[which(aio$budburst == F),])
car::qqPlot(residuals(model))
shapiro.test(residuals(model)) # p < 0.5; residuals not normally distributed

# test for normal distribution for each group individually
library(ggpubr)
ggqqplot(aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),], "ndvi_mean", facet.by = "platform")




## Anova for all platforms and "budburst"
aio %>% 
  filter(budburst == T) %>% 
  filter(platform %in% c("orthomosaic", "planetscope", "sentinel2")) %>% 
  ggplot(aes(x=platform, y=ndvi_mean, group=platform)) +
  geom_boxplot()

anova <- aov(ndvi_mean~platform, data = aio[which(aio$budburst == T & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),])
summary(anova)[[1]] #significantly differ

## Post Hoc Test; compare each group to each other; Tukey procedure; see https://www.beratung-statistik.de/statistik-beratung-infos/r-tutorial/r-varianzanalyse-post-hoc/
# diff - difference of mean values; 
TukeyHSD(aov(ndvi_mean~platform, data = aio[which(aio$budburst == T & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),])) #planetscope - orthomosaic do not differ significantly

## Anova for all platforms and "no budburst"
aio %>% 
  filter(budburst == F) %>% 
  filter(platform %in% c("orthomosaic", "planetscope", "sentinel2")) %>% 
  ggplot(aes(x=platform, y=ndvi_mean, group=platform)) +
  geom_boxplot()

anova <- aov(ndvi_mean~platform, data = aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),])
summary(anova)[[1]] #significantly differ

TukeyHSD(aov(ndvi_mean~platform, data = aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),])) 



