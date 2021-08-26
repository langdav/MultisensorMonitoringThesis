#author: David Langenohl
#last modified: 23.08.2021
#description: derive model performance parameters; MBE (mean bias error), MAD (mean absolute deviation); R²
#NOTE:

rm(list = ls())
library(dplyr)

#load fitted models and manually observed budburst data
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/median_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")


#calculate MBE and MAD
model_fitting_out_mean$MBE <- model_fitting_out_mean$SOS - model_fitting_out_mean$budburst_obervation_doy
model_fitting_out_median$MBE <- model_fitting_out_median$SOS - model_fitting_out_median$budburst_obervation_doy
model_fitting_out_sen1$MBE <- model_fitting_out_sen1$SOS - model_fitting_out_sen1$budburst_obervation_doy

mbe_mad <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  mbe_tmp_mean <- round(sum(model_fitting_out_mean$MBE[model_fitting_out_mean$platform==platform], na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mad_tmp_mean <- round(sum(abs(model_fitting_out_mean$MBE[model_fitting_out_mean$platform==platform]), na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mbe_tmp_median <- round(sum(model_fitting_out_median$MBE[model_fitting_out_median$platform==platform], na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mad_tmp_median <- round(sum(abs(model_fitting_out_median$MBE[model_fitting_out_median$platform==platform]), na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "mean",
                                       mbe = mbe_tmp_mean,
                                       mad = mad_tmp_mean))
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "median",
                                       mbe = mbe_tmp_median,
                                       mad = mad_tmp_median))
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
}
lin_performance <- rbind(lin_performance, data.frame(platform = "sentinel1",
                                                     ndvi_mean_median = NA,
                                                     r_squared = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_sen1))$r.squared,
                                                     slope_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_sen1))$coefficients[8],
                                                     intercept_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_sen1))$coefficients[7]))

plot(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform=="treetalker",])
abline(lm(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform=="treetalker",]))
