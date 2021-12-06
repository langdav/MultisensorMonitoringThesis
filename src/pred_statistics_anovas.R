#author: David Langenohl
#last modified: 15.09.2021
#description: derive model performance; test for differences between estimated and observed budburst dates 
#NOTE:

rm(list = ls())
library(dplyr);library(ggpubr);library(car);library(rstatix)

#load fitted models and manually observed budburst data
# load("out/log_function_models/all_values_fitted_models_output.RData")
# load("out/log_function_models/median_fitted_models_output.RData")
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")
load("out/log_function_models/budburst_fitted_models_output.RData")

#outlier removal (see model_performance_mbe_mad.R)
model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform=="sentinel2" & model_fitting_out_mean$tree_id %in% c("mof_cst_00007",
                                                                                                                      "mof_cst_00026",
                                                                                                                      "mof_cst_00039"))] <- NA #models not fitted right
model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform=="treetalker" & model_fitting_out_mean$tree_id %in% c("mof_cst_00021",
                                                                                                                       "mof_cst_00002",
                                                                                                                       "mof_cst_00012",
                                                                                                                       "mof_cst_00026"))] <- NA #models not fitted right
model_fitting_out_sen1$SOS[model_fitting_out_sen1$tree_id %in% c("mof_cst_00024",
                                                                 "mof_cst_00020",
                                                                 "mof_cst_00001",
                                                                 "mof_cst_00018",
                                                                 "mof_cst_00029",
                                                                 "mof_cst_00030",
                                                                 "mof_cst_00040",
                                                                 "mof_cst_00049")] <- NA #models not fitted right

#add Sentinel-1 data to the other datasets, to be able to compare them with one another
# model_fitting_out_all <- rbind(model_fitting_out_all,model_fitting_out_sen1)
# model_fitting_out_median <- rbind(model_fitting_out_median,model_fitting_out_sen1)
model_fitting_out_mean <- rbind(model_fitting_out_mean,model_fitting_out_sen1)


#check for normality within different groups: platform, estimated - platform, observed
#estimated values; all values
# ggpubr::ggqqplot(model_fitting_out_all,"SOS", facet.by = "platform")
# 
# for(platform in unique(model_fitting_out_all$platform)){
#   pval <- shapiro.test(model_fitting_out_all$SOS[model_fitting_out_all$platform==platform])$p.value
#   print(paste0(platform,": ",pval))
# } 

#estimated values; NDVI means
ggpubr::ggqqplot(model_fitting_out_mean,"SOS", facet.by = "platform")

for(platform in unique(model_fitting_out_mean$platform)){
  pval <- shapiro.test(model_fitting_out_mean$SOS[model_fitting_out_mean$platform==platform])$p.value
  print(paste0(platform,": ",pval))
}

#estimated values; NDVI medians
# ggpubr::ggqqplot(model_fitting_out_median,"SOS", facet.by = "platform")
# 
# for(platform in unique(model_fitting_out_median$platform)){
#   pval <- shapiro.test(model_fitting_out_median$SOS[model_fitting_out_median$platform==platform])$p.value
#   print(paste0(platform,": ",pval))
# } 

#observed values; phase_d
ggpubr::ggqqplot(budburst_model_fitting_out,"SOS", facet.by = "phase")

for(phase in unique(budburst_model_fitting_out$phase)){
  pval <- shapiro.test(budburst_model_fitting_out$SOS[budburst_model_fitting_out$phase==phase])$p.value
  print(paste0(platform,": ",pval))
}

#most of the data is not normally distributed; performing t-tests is not allowed
#Next: test for prerequisites for the performance of an ANOVA
#1. prerequisite: normal distribution; as shown above, this is NOT given; BUT: as the samples are >25, this should not be problematic
#2. prerequisite: variance homogenity
#estimated values; all values
# car::leveneTest(SOS~platform, data = model_fitting_out_all) # p < 0.05; no Variance homogeneity -> Anova can not be used; use Welch-Anova instead

#estimated values; NDVI means
car::leveneTest(SOS~platform, data = model_fitting_out_mean) # p < 0.05; no Variance homogeneity -> Anova can not be used; use Welch-Anova instead

#estimated values; NDVI medians
# car::leveneTest(SOS~platform, data = model_fitting_out_median) # p < 0.05; no Variance homogeneity -> Anova can not be used; use Welch-Anova instead

#observed values
car::leveneTest(SOS~phase, data = budburst_model_fitting_out) # p < 0.05; no Variance homogeneity

#Welch-Anova(s) to test, if estimated values of platforms differ
#estimated values; all values
# model_fitting_out_all %>%
#   group_by() %>% 
#   select(SOS, platform) %>% 
#   rstatix::welch_anova_test(SOS ~ platform)
# p < 0.05 -> groups are significantly different

#estimated values; NDVI means
model_fitting_out_mean %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::welch_anova_test(SOS ~ platform)
# p < 0.05 -> groups are significantly different

#estimated values; NDVI medians
# model_fitting_out_median %>%
#   group_by() %>% 
#   select(SOS, platform) %>% 
#   rstatix::welch_anova_test(SOS ~ platform)
# p < 0.05 -> groups are significantly different

## Tukey post-hoc test; all values
# model_fitting_out_all %>%
#   group_by() %>% 
#   select(SOS, platform) %>% 
#   rstatix::tukey_hsd(SOS ~ platform)

## Tukey post-hoc test; NDVI means
model_fitting_out_mean %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::tukey_hsd(SOS ~ platform)

## Tukey post-hoc test; NDVI median
# model_fitting_out_median %>%
#   group_by() %>% 
#   select(SOS, platform) %>% 
#   rstatix::tukey_hsd(SOS ~ platform)

#Welch-Anova(s) to test, if per-platform estimated values differ from observed values
welch_anova <- NULL
for(platform in unique(model_fitting_out_mean$platform)){

  tmp_df_mean_d <- data.frame(values = c(model_fitting_out_mean$SOS[model_fitting_out_mean$platform==platform],
                                         model_fitting_out_mean$SOS_phase_d[model_fitting_out_mean$platform==platform]),
                              estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,])),
                                                     rep(as.factor("observed"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))))
  
  tmp_df_mean_e <- data.frame(values = c(model_fitting_out_mean$SOS[model_fitting_out_mean$platform==platform],
                                         model_fitting_out_mean$SOS_phase_e[model_fitting_out_mean$platform==platform]),
                              estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,])),
                                                     rep(as.factor("observed"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))))
  
  tmp_df_mean_f <- data.frame(values = c(model_fitting_out_mean$SOS[model_fitting_out_mean$platform==platform],
                                         model_fitting_out_mean$SOS_phase_f[model_fitting_out_mean$platform==platform]),
                              estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,])),
                                                     rep(as.factor("observed"), nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))))
  

  welch_anova <- rbind(welch_anova, data.frame(platform = platform,
                                               mean_median = "mean_d",
                                               p_val = rstatix::welch_anova_test(tmp_df_mean_d, values ~ estimated_observed)$p,
                                               f = as.numeric(rstatix::welch_anova_test(tmp_df_mean_d, values ~ estimated_observed)$statistic)))
                       
  
  welch_anova <- rbind(welch_anova, data.frame(platform = platform,
                                               mean_median = "mean_e",
                                               p_val = rstatix::welch_anova_test(tmp_df_mean_e, values ~ estimated_observed)$p,
                                               f = as.numeric(rstatix::welch_anova_test(tmp_df_mean_e, values ~ estimated_observed)$statistic)))
  
  welch_anova <- rbind(welch_anova, data.frame(platform = platform,
                                               mean_median = "mean_f",
                                               p_val = rstatix::welch_anova_test(tmp_df_mean_f, values ~ estimated_observed)$p,
                                               f = as.numeric(rstatix::welch_anova_test(tmp_df_mean_f, values ~ estimated_observed)$statistic)))
  
}

welch_anova$p_val <- format(round(welch_anova$p_val,6), scientific=F)
