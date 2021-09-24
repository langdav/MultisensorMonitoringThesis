#author: David Langenohl
#last modified: 23.08.2021
#description: derive model performance parameters; MBE (mean bias error), MAD (mean absolute deviation); R²
#NOTE:

rm(list = ls())
library(dplyr)

#load fitted models and manually observed budburst data
load("out/log_function_models/all_values_fitted_models_output.RData")
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/median_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")


#calculate MBE and MAD
model_fitting_out_all$MBE <- model_fitting_out_all$SOS - model_fitting_out_all$budburst_obervation_doy
model_fitting_out_mean$MBE <- model_fitting_out_mean$SOS - model_fitting_out_mean$budburst_obervation_doy
model_fitting_out_median$MBE <- model_fitting_out_median$SOS - model_fitting_out_median$budburst_obervation_doy
model_fitting_out_sen1$MBE <- model_fitting_out_sen1$SOS - model_fitting_out_sen1$budburst_obervation_doy

mbe_mad <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  mbe_tmp_mean <- round(sum(model_fitting_out_mean$MBE[model_fitting_out_mean$platform==platform], na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mad_tmp_mean <- round(sum(abs(model_fitting_out_mean$MBE[model_fitting_out_mean$platform==platform]), na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mbe_tmp_median <- round(sum(model_fitting_out_median$MBE[model_fitting_out_median$platform==platform], na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mad_tmp_median <- round(sum(abs(model_fitting_out_median$MBE[model_fitting_out_median$platform==platform]), na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mbe_tmp_all <- round(sum(model_fitting_out_all$MBE[model_fitting_out_all$platform==platform], na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  mad_tmp_all <- round(sum(abs(model_fitting_out_all$MBE[model_fitting_out_all$platform==platform]), na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  
  
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "mean",
                                       mbe = mbe_tmp_mean,
                                       mad = mad_tmp_mean))
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "median",
                                       mbe = mbe_tmp_median,
                                       mad = mad_tmp_median))
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "all_values",
                                       mbe = mbe_tmp_all,
                                       mad = mad_tmp_all))
}
mbe_tmp_sen1 <- round(sum(model_fitting_out_sen1$MBE, na.rm = T) / nrow(model_fitting_out_sen1),2)
mad_tmp_sen1 <- round(sum(abs(model_fitting_out_sen1$MBE), na.rm = T) / nrow(model_fitting_out_sen1),2)
mbe_mad <- rbind(mbe_mad, data.frame(platform = "Sentinel-1",
                                     ndvi_mean_median = NA,
                                     mbe = mbe_tmp_sen1,
                                     mad = mad_tmp_sen1))



#perform a simple linear regression between estimated and observed budburst dates and get the rsquared of the model
lin_performance <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       ndvi_mean_median = "mean",
                                                       r_squared = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,
                                                       slope_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[8],
                                                       intercept_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[7]))
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       ndvi_mean_median = "median",
                                                       r_squared = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$r.squared,
                                                       slope_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[8],
                                                       intercept_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[7]))
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       ndvi_mean_median = "all_values",
                                                       r_squared = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$r.squared,
                                                       slope_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[8],
                                                       intercept_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[7]))
}
lin_performance <- rbind(lin_performance, data.frame(platform = "sentinel1",
                                                     ndvi_mean_median = NA,
                                                     r_squared = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_sen1))$r.squared,
                                                     slope_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_sen1))$coefficients[8],
                                                     intercept_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_sen1))$coefficients[7]))

lin_performance$r_squared <- format(round(lin_performance$r_squared,3), scientific=F)

# manually check for distribution of points in the fitted R² model; outliers?; few points?
#mean per-tree NDVI values
ggplot(model_fitting_out_mean, 
       aes(x=SOS, y=budburst_obervation_doy)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year")

#median per-tree NDVI values
ggplot(model_fitting_out_median, 
       aes(x=SOS, y=budburst_obervation_doy)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year")

#all per-tree NDVI values
ggplot(model_fitting_out_all, 
       aes(x=SOS, y=budburst_obervation_doy)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year")

#sentinel-1
ggplot(model_fitting_out_sen1, 
       aes(x=SOS, y=budburst_obervation_doy)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  theme_light() +
  xlab("Day of Year")

#find earliest budburst, latest budburst and total number of detected budbursts per platform
platform <- "orthomosaic"

min(model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform == platform)], na.rm = T)
max(model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform == platform)], na.rm = T)
length(which(is.na(model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform == platform)])==F))

min(model_fitting_out_median$SOS[which(model_fitting_out_median$platform == platform)], na.rm = T)
max(model_fitting_out_median$SOS[which(model_fitting_out_median$platform == platform)], na.rm = T)
length(which(is.na(model_fitting_out_median$SOS[which(model_fitting_out_median$platform == platform)])==F))

min(model_fitting_out_all$SOS[which(model_fitting_out_all$platform == platform)], na.rm = T)
max(model_fitting_out_all$SOS[which(model_fitting_out_all$platform == platform)], na.rm = T)
length(which(is.na(model_fitting_out_all$SOS[which(model_fitting_out_all$platform == platform)])==F))

# Sentinel-1
min(model_fitting_out_sen1$SOS, na.rm = T)
max(model_fitting_out_sen1$SOS, na.rm = T)
length(which(is.na(model_fitting_out_sen1$SOS)==F))






which(is.na(model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform == "treetalker")])==T)
model_fitting_out_mean[which(model_fitting_out_mean$platform == "treetalker"),]
