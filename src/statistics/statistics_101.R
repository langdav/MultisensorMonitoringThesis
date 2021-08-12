rm(list = ls())
library(plyr);library(dplyr);library(ggplot2)

## load sample data
tt_data <- readRDS("data/TreeTalker/TT_spectraldata.RDS")
tt_data$AS7263_610 <- as.numeric(tt_data$AS7263_610)
tt_data <- tt_data %>% filter(tt_data[,c(5:16)] > 5) ## remove entries < 5, as they are not needed for the analysis
#head(tt_data)

#########################################################
## check for outliers ##
#######################

# boxplot for visual detection of outliers
boxplot(tt_data$AS7263_610)
boxplot(log(tt_data$AS7263_610)) #log-transformed

###################
# which values are outliers (based on Inter Quantile Range: Q25 - 1.5 * IQR  <-> Q75 + 1.5*IQR; IQR = Q75 - Q25)
quantile(tt_data$AS7263_610, 0.25) #first quartile
quantile(tt_data$AS7263_610, 0.75) #third quartile
IQR(tt_data$AS7263_610) #innerquartile range

boxplot.stats(tt_data$AS7263_610)$out # values wich lay beyond the interquartile  range
out_ind <- which(tt_data$AS7263_610 %in% c(boxplot.stats(tt_data$AS7263_610)$out)) # locate outlier rows
AS7263_610_clean <- tt_data$AS7263_610[-out_ind] #remove outlier rows

# same for log-transformed values
out_ind <- which(log(tt_data$AS7263_610) %in% c(boxplot.stats(log(tt_data$AS7263_610))$out)) # locate outlier rows
AS7263_610_log_clean <- log(tt_data$AS7263_610[-out_ind]) #remove outlier rows

####################
# alternative method: find outliers according to 2.5 and 97.5 percentiles
quantile(tt_data$AS7263_610, 0.025) #lower quantile; all values < this value are outliers
quantile(tt_data$AS7263_610, 0.975) #upper quantile; all values > this value are outliers

out_ind <- which(tt_data$AS7263_610 < quantile(tt_data$AS7263_610, 0.025) | tt_data$AS7263_610 > quantile(tt_data$AS7263_610, 0.975)) # locate outlier rows
AS7263_610_clean <- tt_data$AS7263_610[-out_ind] #remove outlier rows

# same for log-transformed values
out_ind <- which(log(tt_data$AS7263_610) < quantile(log(tt_data$AS7263_610), 0.025) | log(tt_data$AS7263_610) > quantile(log(tt_data$AS7263_610), 0.975)) # locate outlier rows
AS7263_610_log_clean <- log(tt_data$AS7263_610[-out_ind]) #remove outlier rows

####################
# alternative method: identify outliers using the rstatix package
# example data: all-in-one data frame with ndvi values from all platforms
load("out/all_in_one/all_platforms_daily_ndvi_values.RData")

# we are looking for outliers within each individual group "orthomosaic", "planetscope", "sentinel2"
library(rstatix)
aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),] %>% 
  group_by(platform) %>%
  identify_outliers(ndvi_mean)

# we are looking for outliers within each individual group "orthomosaic", "planetscope", "sentinel2" and for each individual tree
aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),] %>% 
  group_by(platform) %>%
  group_by(tree_id) %>% 
  identify_outliers(ndvi_mean)

#########################################################
## check for normality ##
########################
  
################
## graphically
################

## histogram of single spectrum with normal distribution curve overlay 
#par(mfrow = c(1,1))
data <- tt_data$AS7263_610
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F, main = paste(colnames(data)))
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

## histogram of LOG TRANSFORMED single spectrum with normal distribution curve overlay
data <- log(tt_data$AS7263_610)
std <- sd(data,na.rm=T) #standard deviation of data; for normal distribution
m <- mean(data,na.rm=T) #mean of data; for normal distribution
hist(data, freq=F, main = paste(colnames(data)))
curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red")

## multiple histograms with normal distribution curve overlay
data <- log(tt_data[,5:16])
std <- lapply(data, sd, na.rm = T)
m <- lapply(data, mean, na.rm = T)

par(mfrow=c(4,3))
for(i in 1:12) {
  hist(data[,i], freq=F, main = paste(colnames(data)[i]))
  curve(dnorm(x,mean=m[[i]],sd=std[[i]]), add=TRUE,col="red")
}

## QQ-Plot
qqPlot(tt_data$AS7263_610) #not normally distributed
qqPlot(log(tt_data$AS7263_610)) #not normally distributed
qqPlot(sqrt(tt_data$AS7263_610)) #not normally distributed

## multiple QQ-Plots
library(car) #for qqPlot-function

data <- log(tt_data[,5:16])
par(mfrow=c(4,3))
for(i in 1:12) {
  qqPlot(data[,i])
}

####################
# alternative method: multiple QQ-Plots with package ggpubr
# example data: all-in-one data frame with ndvi values from all platforms
load("out/all_in_one/all_platforms_daily_ndvi_values.RData")

# we are looking for normality within each individual group "orthomosaic", "planetscope", "sentinel2"
library(ggpubr)
ggqqplot(aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),], "ndvi_mean", facet.by = "platform")

####################
## statistically
###################

## shapiro-Wilk test for normal distribution; only applicable for data < 5000
# a p-value > 0.05 means that the data is normally distributed (as H0: data IS normally distributed)
shapiro.test(tt_data$AS7263_610) # dataset too big (> 5000); another test is needed (e.g. a graphical test)

##################################################################
## check for variance homogenity (Anova prerequisite) ##
#######################################################
# example data: all-in-one data frame with ndvi values from all platforms
load("out/all_in_one/all_platforms_daily_ndvi_values.RData")

library(car)
leveneTest(ndvi_mean~platform, data = aio[which(aio$budburst == F & aio$platform %in% c("orthomosaic", "planetscope", "sentinel2")),]) 
# p < 0.05; no Variance homogeneity -> Anova can not be used; use Welch-Anova instead

##################################################################
## check for (and remove) outliers ##
####################################

###################################
## Boxplot and Quantiles/Quartiles
boxplot(tt_data$AS7263_610)
boxplot(log(tt_data$AS7263_610)) #log transformed

## which values are the outliers (based on Inter Quantile Range: Q25 - 1.5 * IQR  <-> Q75 + 1.5*IQR; IQR = Q75 - Q25)
quantile(tt_data$AS7263_610, 0.25) #first quartile
quantile(tt_data$AS7263_610, 0.75) #third quartile
IQR(tt_data$AS7263_610) #innerquartile range

boxplot.stats(tt_data$AS7263_610)$out # values which do not lay within the interquartile range
out_ind <- which(tt_data$AS7263_610 %in% c(boxplot.stats(tt_data$AS7263_610)$out)) 
out_ind #rows containing the outliers
outliers_removed <- tt_data$AS7263_610[-out_ind] #remove rows that are containing the outliers

## alternative method: find outliers according to 2.5 and 97.5 percentiles
quantile(tt_data$AS7263_610, 0.025) #lower quantile; all values < this value are outliers
quantile(tt_data$AS7263_610, 0.975) #upper quantile; all values > this value are outliers

out_ind <- which(tt_data$AS7263_610 < quantile(tt_data$AS7263_610, 0.025) | tt_data$AS7263_610 > quantile(tt_data$AS7263_610, 0.975))
out_ind #rows containing the outliers
outliers_removed <- tt_data$AS7263_610[-out_ind] #remove rows that are containing the outliers

## after removing outliers: recheck for normality
# outliers removed with ITR
car::qqPlot(outliers_removed) #nope
car::qqPlot(log(outliers_removed)) #log-transformed; not really

std <- sd(outliers_removed,na.rm=T) 
m <- mean(outliers_removed,na.rm=T)
hist(outliers_removed, freq=F, main = paste(colnames(outliers_removed)));curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red") #nope

std <- sd(log(outliers_removed),na.rm=T) 
m <- mean(log(outliers_removed),na.rm=T)
hist(log(outliers_removed), freq=F, main = paste(colnames(log(outliers_removed))));curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red") #log-transformed; not really

shapiro.test(outliers_removed) #dataset still to big

# outliers removed with 2.5 and 97.5 percentiles
car::qqPlot(outliers_removed) #nope
car::qqPlot(log(outliers_removed)) #log-transformed; not really

std <- sd(outliers_removed,na.rm=T) 
m <- mean(outliers_removed,na.rm=T)
hist(outliers_removed, freq=F, main = paste(colnames(outliers_removed)));curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red") #nope

std <- sd(log(outliers_removed),na.rm=T) 
m <- mean(log(outliers_removed),na.rm=T)
hist(log(outliers_removed), freq=F, main = paste(colnames(log(outliers_removed))));curve(dnorm(x,mean=m,sd=std), add=TRUE,col="red") #log-transformed; not really

##################################################################
## testing for differences between groups ##
###########################################

# test data: Sentinel-2 NDVI values
sen2_ndvi <- read.csv("out/sentinel2/ndvi_long_format_phenoclasses_sentinel2.csv")
sen2_ndvi$date <- as.Date(sen2_ndvi$date)

# test data: Planetscope NDVI values
plan_ndvi <- read.csv("out/planetscope/ndvi_long_format_phenoclasses_planetscape.csv")
plan_ndvi$date <- as.Date(plan_ndvi$date)
names(plan_ndvi) <- c("tree_id","date","ndvi","budburst","budburst_perc")

#############################
## t-tests; comparing if two means are significantly different; H0 = "they are equal" so p<0.05 means, they are different
# requirements:
# Two independent samples/groups
# metric scaled y-variable
# normally distributed y-variable within the groups
# Homogeneous (almost equal) variances of the y-variables of the groups (Levene test)
# In the case of unequal variances, the so-called Welch test or Welch-t test is calculated

t.test(x=sen2_ndvi$ndvi[which(sen2_ndvi$budburst == T)], y=sen2_ndvi$ndvi[which(sen2_ndvi$budburst == F)]) # p < 0.05; different
# note: this test was just performed for demonstration purposes; those variables are NOT normally distributed, so the t-test is not allowed here

#############################
## Wilcoxon-Test; comparing if two means are significantly different; y-values do NOT need to be normally distributed
# requirements:
# two independent samples/groups – three or more groups via the Kruskal-Wallis test
# ordinal or metric scaled y-variable
# normally distributed y-variable within the groups NOT necessary

wilcox.test(ndvi~budburst, data=sen2_ndvi,exact=F) #NDVI grouped by budburst ("budburst" and "no budburst");p < 0.05 <- different
z <- qnorm(wilcox.test(ndvi~budburst, data=sen2_ndvi,exact=F)$p.value)
r <- abs(z/sqrt(nrow(sen2_ndvi)))
r # > 0.3 aber < 0.5 -> mittlere Effektstärke des budburst auf ndvi

#############################
## Spearman correlation coefficient; ordinal scaled variable (budburst percent) ~ metric scaled variable (NDVI); tests for correlation
# requirements:
# TWO ordinal scaled variables or one metric scaled and one ordinal scaled variable

cor.test(sen2_ndvi$budburst_perc, sen2_ndvi$ndvi, method = "spearman") # p < 0.05 -> correlation; rho > 0.3 und < 0.5 -> mittlere Effektstärke
cor.test(sen2_ndvi$budburst_perc, sen2_ndvi$ndvi, method = "spearman", alternative = "greater") # if alternative = "greater": p < 0.05 -> positive correlation
cor.test(sen2_ndvi$budburst_perc, sen2_ndvi$ndvi, method = "spearman", alternative = "less") # if alternative = "greater": p < 0.05 -> negative correlation

#############################
## ANOVA
# requirements: 
# More THAN TWO independent samples/groups
# metric scaled y-variable
# Normally distributed error terms within the groups
# Homogeneous (almost equal) variances of the y-variables of the groups (descriptive or Levene test)

## test for normal distribution of the residuals (error terms)
# linear model; the residuals are then taken from it
model1 <- lm(ndvi~budburst_perc, data = plan_ndvi)
model2 <- lm(ndvi~budburst_perc, data = sen2_ndvi)

# qqplots
car::qqPlot(rstandard(model1)) # NOT normally distributed residuals
car::qqPlot(rstandard(model2)) #somewhat normally distributed residuals

## Anova
anova_plan <- aov(plan_ndvi$ndvi~plan_ndvi$budburst_perc)
summary(anova_plan) # p < 0.05 -> means differ

anova_sen2 <- aov(sen2_ndvi$ndvi~sen2_ndvi$budburst_perc)
summary(anova_sen2) # p < 0.05 -> means differ

# pairwise t-test (every group against every other group)
pairwise.t.test(plan_ndvi$ndvi,plan_ndvi$budburst_perc, p.adjust="bonferroni")
pairwise.t.test(sen2_ndvi$ndvi,sen2_ndvi$budburst_perc, p.adjust="bonferroni")
