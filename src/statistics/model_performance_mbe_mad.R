#author: David Langenohl
#last modified: 23.08.2021
#description: derive model performance parameters; MBE (mean bias error), MAD (mean absolute deviation); R²
#NOTE:

rm(list = ls())
#library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);
library(dplyr)

#load fitted models and manually observed budburst data
load("out/log_function_models/fitted_models_output.RData")

#calculate MBE and MAD
model_fitting_out$MBE <- model_fitting_out$SOS - model_fitting_out$budburst_obervation_doy

mbe_mad <- NULL
for(platform in unique(model_fitting_out$platform)){
  mbe_tmp <- round(sum(model_fitting_out$MBE[model_fitting_out$platform==platform], na.rm = T) / nrow(model_fitting_out[model_fitting_out$platform==platform,]),2)
  mad_tmp <- round(sum(abs(model_fitting_out$MBE[model_fitting_out$platform==platform]), na.rm = T) / nrow(model_fitting_out[model_fitting_out$platform==platform,]),2)
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                               mbe = mbe_tmp,
                               mad = mad_tmp))
}

#perform a simple linear regression between estimated and observed budburst dates and get the rsquared of the model
lin_performance <- NULL
for(platform in unique(model_fitting_out$platform)){
  rsquared_tmp <- summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out[model_fitting_out$platform==platform,]))$r.squared #R² ~ 34%; p = signif slope; NOT singif intercept
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       r_squared = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out[model_fitting_out$platform==platform,]))$r.squared,
                                                       slope_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out[model_fitting_out$platform==platform,]))$coefficients[8],
                                                       intercept_signif = summary(lm(SOS ~ budburst_obervation_doy,model_fitting_out[model_fitting_out$platform==platform,]))$coefficients[7]))
}

plot(SOS ~ budburst_obervation_doy,model_fitting_out[model_fitting_out$platform=="treetalker",])
abline(lm(SOS ~ budburst_obervation_doy,model_fitting_out[model_fitting_out$platform=="treetalker",]))
