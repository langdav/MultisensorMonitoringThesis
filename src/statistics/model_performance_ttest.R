#author: David Langenohl
#last modified: 15.09.2021
#description: derive model performance; test for differences between estimated and observed budburst dates 
#NOTE:

rm(list = ls())
library(dplyr);library(ggpubr);library(car);library(rstatix)

#load fitted models and manually observed budburst data
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/median_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")
model_fitting_out_sen1$platform <- "sentinel1"


#add Sentinel-1 data to the other datasets, to be able to compare them with one another
model_fitting_out_mean <- rbind(model_fitting_out_mean,model_fitting_out_sen1)
model_fitting_out_median <- rbind(model_fitting_out_median,model_fitting_out_sen1)

#check for normality within different groups: platform, estimated - platform, observed
#estimated values; NDVI means
ggpubr::ggqqplot(model_fitting_out_mean,"SOS", facet.by = "platform")

for(platform in unique(model_fitting_out_mean$platform)){
  pval <- shapiro.test(model_fitting_out_mean$SOS[model_fitting_out_mean$platform==platform])$p.value
  print(paste0(platform,": ",pval))
} #nothing normally distributed

#estimated values; NDVI medians
ggpubr::ggqqplot(model_fitting_out_median,"SOS", facet.by = "platform")

for(platform in unique(model_fitting_out_median$platform)){
  pval <- shapiro.test(model_fitting_out_median$SOS[model_fitting_out_median$platform==platform])$p.value
  print(paste0(platform,": ",pval))
} #Planetscope and Orthomosaic normally distributed; Treetalker and Sentinel-2 not normally distributed

#observed values
ggpubr::ggqqplot(model_fitting_out_mean,"budburst_obervation_doy", facet.by = "platform")

for(platform in unique(model_fitting_out_mean$platform)){
  pval <- shapiro.test(model_fitting_out_mean$budburst_obervation_doy[model_fitting_out_mean$platform==platform])$p.value
  print(paste0(platform,": ",pval))
}

#NONE of the data is normally distributed; t-test performance is not allowed
#Next: test for prerequisites for the performance of an ANOVA
#1. prerequisite: normal distribution; as shown above, this is NOT given; BUT: as the samples are >25, this should not be problematic
#2. prerequisite: variance homogenity
#estimated values; NDVI means
car::leveneTest(SOS~platform, data = model_fitting_out_mean) # p < 0.05; no Variance homogeneity -> Anova can not be used; use Welch-Anova instead

#estimated values; NDVI medians
car::leveneTest(SOS~platform, data = model_fitting_out_median) # p < 0.05; no Variance homogeneity -> Anova can not be used; use Welch-Anova instead

#observed values
car::leveneTest(budburst_obervation_doy~platform, data = model_fitting_out_mean) # p > 0.05; Variance homogeneity -> Anova can be used

#Welch-Anova(s) to test, if estimated values of platforms differ
#estimated values; NDVI means
model_fitting_out_mean %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::welch_anova_test(SOS ~ platform)
# p < 0.05 -> groups are significantly different

#estimated values; NDVI medians
model_fitting_out_median %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::welch_anova_test(SOS ~ platform)
# p < 0.05 -> groups are significantly different

## Tukey post-hoc test; NDVI means
model_fitting_out_mean %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::tukey_hsd(SOS ~ platform)

## Tukey post-hoc test; NDVI median
model_fitting_out_median %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::tukey_hsd(SOS ~ platform)

#Welch-Anova(s) to test, if per-platform estimated values differ from observed values
welch_anova <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  tmp_df_mean <- data.frame(values = c(model_fitting_out_mean$SOS[model_fitting_out_mean$platform==platform],
                               model_fitting_out_mean$budburst_obervation_doy[model_fitting_out_mean$platform==platform]),
                    estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,])),
                                           rep(as.factor("observed"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))))
  
  tmp_df_median <- data.frame(values = c(model_fitting_out_median$SOS[model_fitting_out_median$platform==platform],
                                       model_fitting_out_median$budburst_obervation_doy[model_fitting_out_median$platform==platform]),
                            estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,])),
                                                   rep(as.factor("observed"), nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]))))
  
  welch_anova <- rbind(welch_anova, data.frame(platform = platform,
                                               mean_median = "mean",
                                       p_val = rstatix::welch_anova_test(tmp_df_mean, values ~ estimated_observed)$p))
  
  welch_anova <- rbind(welch_anova, data.frame(platform = platform,
                                               mean_median = "median",
                                               p_val = rstatix::welch_anova_test(tmp_df_median, values ~ estimated_observed)$p))
}
