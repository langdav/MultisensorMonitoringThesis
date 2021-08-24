#author: David Langenohl
#last modified: 23.08.2021
#description: derive model performance parameters; MBE (mean bias error), MAD (mean absolute deviation); R²
#NOTE:

rm(list = ls())
library(dplyr)

#load fitted models and manually observed budburst data
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/median_fitted_models_output.RData")

#calculate MBE and MAD
model_fitting_out_mean$MBE <- model_fitting_out_mean$SOS - model_fitting_out_mean$budburst_obervation_doy
model_fitting_out_median$MBE <- model_fitting_out_median$SOS - model_fitting_out_median$budburst_obervation_doy

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

#perform a simple linear regression between estimated and observed budburst dates and get the rsquared of the model
lin_performance <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  rsquared_tmp_mean <- summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared #R² ~ 34%; p = signif slope; NOT singif intercept
  rsquared_tmp_median <- summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$r.squared #R² ~ 34%; p = signif slope; NOT singif intercept
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


plot(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform=="treetalker",])
abline(lm(SOS ~ budburst_obervation_doy,model_fitting_out_mean[model_fitting_out_mean$platform=="treetalker",]))
