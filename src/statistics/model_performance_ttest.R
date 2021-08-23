#author: David Langenoh
#last modified: 23.08.2021
#description: derive model performance; test for differences between estimated and observed budburst dates 
#NOTE:

rm(list = ls())
library(dplyr);library(ggpubr);library(car);library(rstatix)

#load fitted models and manually observed budburst data
load("out/log_function_models/fitted_models_output.RData")

#check for normality within different groups: platform, estimated - platform, observed
#estimated values
ggpubr::ggqqplot(model_fitting_out,"SOS", facet.by = "platform")

for(platform in unique(model_fitting_out$platform)){
  pval <- shapiro.test(model_fitting_out$SOS[model_fitting_out$platform==platform])$p.value
  print(paste0(platform,": ",pval))
}

#observed values
ggpubr::ggqqplot(model_fitting_out,"budburst_obervation_doy", facet.by = "platform")

for(platform in unique(model_fitting_out$platform)){
  pval <- shapiro.test(model_fitting_out$budburst_obervation_doy[model_fitting_out$platform==platform])$p.value
  print(paste0(platform,": ",pval))
}

  
#NONE of the data is normally distributed; t-test performance is not allowed
#Next: test for prerequisites for the performance of an ANOVA
#1. prerequisite: normal distribution; as shown above, this is NOT given; BUT: as the samples are >25, this should not be problematics
#2. prerequisite: variance homogenity
#estimated values
car::leveneTest(SOS~platform, data = model_fitting_out) # p < 0.05; no Variance homogeneity -> Anova can not be used; use Welch-Anova instead

#observed values
car::leveneTest(budburst_obervation_doy~platform, data = model_fitting_out) # p > 0.05; Variance homogeneity -> Anova can be used

#Welch-Anova(s) to test, if estimated values of platforms differ
#estimated values
model_fitting_out %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::welch_anova_test(SOS ~ platform)
# p < 0.05 -> groups are significantly different

## Tukey post-hoc test
model_fitting_out %>%
  group_by() %>% 
  select(SOS, platform) %>% 
  rstatix::tukey_hsd(SOS ~ platform)
# orthomosaic - planetscope; not signif.
# orthomosaic - sentinel2; signif.
# orthomosaic - treetalker; signif.
# planetscope - sentinel2; signif.
# planetscope - treetalker; signif.
# sentinel2 - treetalker; not signif.

#Welch-Anova(s) to test, if per-platform estimated values differ from observed values
welch_anova <- NULL
for(platform in unique(model_fitting_out$platform)){
  tmp_df <- data.frame(values = c(model_fitting_out$SOS[model_fitting_out$platform==platform],
                               model_fitting_out$budburst_obervation_doy[model_fitting_out$platform==platform]),
                    estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out[model_fitting_out$platform==platform,])),
                                           rep(as.factor("observed"), nrow(model_fitting_out[model_fitting_out$platform==platform,]))))
  
  welch_anova <- rbind(welch_anova, data.frame(platform = platform,
                                       p_val = rstatix::welch_anova_test(tmp_df, values ~ estimated_observed)$p))
}