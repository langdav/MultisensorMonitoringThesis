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

##################################################################
## visual comparisons ######
###########################
## NDVI~Budburst; Sentinel-2 and Planetscope; all NDVI values from all trees
# sentinel2 and planetscope in one dataframe for plotting; adding an additional "platform" coloumn
sen2_ndvi$platform <- rep("sentinel2", nrow(sen2_ndvi))
plan_ndvi$platform <- rep("planetscope", nrow(plan_ndvi))
platforms_tmp <- rbind(sen2_ndvi, plan_ndvi)

# reduce to dates that are included in both, sentinel2 and planetscope, data; hence only three days (thanks to sentinel2)
platforms_tmp <- platforms_tmp %>% filter(date %in% as.Date(c("2021-03-20","2021-03-30","2021-05-09")))

# compare both of them
platforms_tmp %>% 
  ggplot(aes(x=as.factor(budburst), y = ndvi, color = platform)) +
  geom_boxplot()

##################################################################
## statistic comparisons ########################################
################################################################

## shapiro-Wilk test for normal distribution of NDVI-values
# Sentinel-2 -> outcome: NDVI normally distributed for budburst = T, but not for budburst = F
shapiro.test(sen2_ndvi$ndvi[which(sen2_ndvi$budburst == T)]) # p > 0.05 -> normalverteilt
shapiro.test(sen2_ndvi$ndvi[which(sen2_ndvi$budburst == F)]) # p < 0.05 -> nicht normalverteilt

# Planetscope -> outcome: NDVI not normally distributed; not for budburst = T, nor for budburst = F
shapiro.test(plan_ndvi$ndvi[which(plan_ndvi$budburst == T)]) # too large; > 5000 points; check graphically

################################
## check for normal distribution graphically (histogram & qqplot)
# Sentinel-2
hist(sen2_ndvi$ndvi[which(sen2_ndvi$budburst == F)]) ## histogram;  not normal distributed
test <- scale(sen2_ndvi$ndvi[which(sen2_ndvi$budburst == F)]);qqnorm(test);qqline(test) #Q-Q-Plot; z-standardisation;clearly not normal distributed!

# Planetscope
hist(plan_ndvi$ndvi[which(plan_ndvi$budburst == T)]) ## histogram; clearly not normal distributed!
test <- scale(plan_ndvi$ndvi[which(plan_ndvi$budburst == T)]);qqnorm(test);qqline(test) #Q-Q-Plot; z-standardisation;clearly not normal distributed!

#############################
## t-tests; comparing if two means are significantly different; H0 = "they are equal" so p<0.05 means, they are different
# Voraussetungen:
# wei voneinander unabhängige Stichproben/Gruppen
# metrisch skalierte y-Variable
# normalverteilte y-Variable innerhalb der Gruppen
# Homogene (nahezu gleiche) Varianzen der y-Variablen der Gruppen (Levene-Test)
# Bei ungleichen Varianzen wird der sog. Welch-Test bzw. Welch-t-Test gerechnet
####### kann ich so an sich nicht anwenden, da die y-Werte innerhalb der Gruppen NICHT (alle) normalverteilt sind; daher: Wilcoxon-Test

# Planetscope; compare groups "budburst true" and "budburst false"
t.test(x=plan_ndvi$ndvi[which(plan_ndvi$budburst == T)], y=plan_ndvi$ndvi[which(plan_ndvi$budburst == F)]) #different

# Sentinel-2; compare groups "budburst true" and "budburst false"
t.test(x=sen2_ndvi$ndvi[which(sen2_ndvi$budburst == T)], y=sen2_ndvi$ndvi[which(sen2_ndvi$budburst == F)]) #different

# Planetscope vs Sentinel-2; budburst = T
t.test(x=plan_ndvi$ndvi[which(plan_ndvi$budburst == T)], y=sen2_ndvi$ndvi[which(sen2_ndvi$budburst == T)]) #different

# Planetscope vs Sentinel-2; budburst = F
t.test(x=plan_ndvi$ndvi[which(plan_ndvi$budburst == F)], y=sen2_ndvi$ndvi[which(sen2_ndvi$budburst == F)]) #different

#############################
## Wilcoxon-Test; comparing if two means are significantly different; y-values do NOT need to be normally distributed
# Voraussetzungen
# zwei voneinander unabhängige Stichproben/Gruppen – drei oder mehr Gruppen über den Kruskal-Wallis-Test
# ordinal oder metrisch skalierte y-Variable
# normalverteilte y-Variable innerhalb der Gruppen nicht nötig

# plan_ndvi$budburst_num <- rep(0, nrow(plan_ndvi))
# plan_ndvi$budburst_num[which(plan_ndvi$budburst == T)] <- 1
wilcox.test(ndvi~budburst, data=plan_ndvi,exact=F) #p < 0.05 <- different
z <- qnorm(wilcox.test(ndvi~budburst, data=plan_ndvi,exact=F)$p.value)
r <- abs(z/sqrt(nrow(plan_ndvi)))
r # ?

wilcox.test(ndvi~budburst, data=sen2_ndvi,exact=F) #p < 0.05 <- different
z <- qnorm(wilcox.test(ndvi~budburst, data=sen2_ndvi,exact=F)$p.value)
r <- abs(z/sqrt(nrow(sen2_ndvi)))
r # > 0.3 aber < 0.5 -> mittlere Effektstärke des budburst auf ndvi


#############################
## Spearman correlation coefficient; ordinal scaled variable (budburst percent) ~ metric scaled variable (NDVI); tests for correlation
cor.test(plan_ndvi$budburst_perc, plan_ndvi$ndvi, method = "spearman") # p < 0.05 -> correlation; rho > 0.5 -> große Effektstärke
cor.test(plan_ndvi$budburst_perc, plan_ndvi$ndvi, method = "spearman", alternative = "greater") # p < 0.05 -> positive correlation
cor.test(plan_ndvi$budburst_perc, plan_ndvi$ndvi, method = "spearman", alternative = "less") # p = 1 -> not a negative correlation

cor.test(sen2_ndvi$budburst_perc, sen2_ndvi$ndvi, method = "spearman") # p < 0.05 -> correlation; rho > 0.3 und < 0.5 -> mittlere Effektstärke
cor.test(sen2_ndvi$budburst_perc, sen2_ndvi$ndvi, method = "spearman", alternative = "greater") # p < 0.05 -> positive correlation
cor.test(sen2_ndvi$budburst_perc, sen2_ndvi$ndvi, method = "spearman", alternative = "less") # p = 1 -> not a negative correlation

#############################
## ANOVA
# Voraussetzungen: 
# mehr als zwei voneinander unabhängige Stichproben/Gruppen
# metrisch skalierte y-Variable
# normalverteilte Fehlerterme innerhalb der Gruppen
# Homogene (nahezu gleiche) Varianzen der y-Variablen der Gruppen (deskriptiv oder Levene-Test)

## testen auf Normalverteilung der Residuen (Fehlerterme)
# lineares modell; davon werden dann die Residuen genommen
model1 <- lm(ndvi~budburst_perc, data = plan_ndvi)
model2 <- lm(ndvi~budburst_perc, data = sen2_ndvi)

# qqplots
qqnorm(rstandard(model1));qqline(rstandard(model1)) # keine Ahnung
qqnorm(rstandard(model2));qqline(rstandard(model2)) # keine Ahnung

## Anova
anova_plan <- aov(plan_ndvi$ndvi~plan_ndvi$budburst_perc)
summary(anova_plan) # p < 0.05 -> means differ

anova_sen2 <- aov(sen2_ndvi$ndvi~sen2_ndvi$budburst_perc)
summary(anova_sen2) # p < 0.05 -> means differ

# pairwise t-test (every group against every other group)
pairwise.t.test(plan_ndvi$ndvi,plan_ndvi$budburst_perc, p.adjust="bonferroni")
pairwise.t.test(sen2_ndvi$ndvi,sen2_ndvi$budburst_perc, p.adjust="bonferroni")

################################
## testing for differences between the NDVI of individual trees within the different products
# Planetscape
anova_plan <- aov(plan_ndvi$ndvi~plan_ndvi$tree_id)
summary(anova_plan) # p < 0.05 -> means differ

test <- pairwise.t.test(plan_ndvi$ndvi,plan_ndvi$tree_id, p.adjust="bonferroni")
View(test$p.value)
