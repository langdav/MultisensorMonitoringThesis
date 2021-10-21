#author: David Langenohl
#last modified: 23.08.2021
#description: derive model performance parameters; MBE (mean bias error), MAD (mean absolute deviation); R²
#NOTE:

rm(list = ls())
library(dplyr);library(ggplot2)

#load fitted models and manually observed budburst data
load("out/log_function_models/all_values_fitted_models_output.RData")
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/median_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")
load("out/log_function_models/budburst_fitted_models_output.RData")

#calculate MBE and MAD
model_fitting_out_all$MBE_d <- model_fitting_out_all$SOS - model_fitting_out_all$SOS_phase_d
model_fitting_out_all$MBE_e <- model_fitting_out_all$SOS - model_fitting_out_all$SOS_phase_e
model_fitting_out_all$MBE_f <- model_fitting_out_all$SOS - model_fitting_out_all$SOS_phase_f
model_fitting_out_mean$MBE_d <- model_fitting_out_mean$SOS - model_fitting_out_mean$SOS_phase_d
model_fitting_out_mean$MBE_e <- model_fitting_out_mean$SOS - model_fitting_out_mean$SOS_phase_e
model_fitting_out_mean$MBE_f <- model_fitting_out_mean$SOS - model_fitting_out_mean$SOS_phase_f
model_fitting_out_median$MBE_d <- model_fitting_out_median$SOS - model_fitting_out_median$SOS_phase_d
model_fitting_out_median$MBE_e <- model_fitting_out_median$SOS - model_fitting_out_median$SOS_phase_e
model_fitting_out_median$MBE_f <- model_fitting_out_median$SOS - model_fitting_out_median$SOS_phase_f
model_fitting_out_sen1$MBE_d <- model_fitting_out_sen1$SOS - model_fitting_out_sen1$SOS_phase_d
model_fitting_out_sen1$MBE_e <- model_fitting_out_sen1$SOS - model_fitting_out_sen1$SOS_phase_e
model_fitting_out_sen1$MBE_f <- model_fitting_out_sen1$SOS - model_fitting_out_sen1$SOS_phase_f

mbe_mad <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  mbe_d_tmp_mean <- round(sum(model_fitting_out_mean$MBE_d[model_fitting_out_mean$platform==platform], na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mbe_e_tmp_mean <- round(sum(model_fitting_out_mean$MBE_e[model_fitting_out_mean$platform==platform], na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mbe_f_tmp_mean <- round(sum(model_fitting_out_mean$MBE_f[model_fitting_out_mean$platform==platform], na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mad_d_tmp_mean <- round(sum(abs(model_fitting_out_mean$MBE_d[model_fitting_out_mean$platform==platform]), na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mad_e_tmp_mean <- round(sum(abs(model_fitting_out_mean$MBE_e[model_fitting_out_mean$platform==platform]), na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  mad_f_tmp_mean <- round(sum(abs(model_fitting_out_mean$MBE_f[model_fitting_out_mean$platform==platform]), na.rm = T) / nrow(model_fitting_out_mean[model_fitting_out_mean$platform==platform,]),2)
  
  mbe_d_tmp_median <- round(sum(model_fitting_out_median$MBE_d[model_fitting_out_median$platform==platform], na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mbe_e_tmp_median <- round(sum(model_fitting_out_median$MBE_e[model_fitting_out_median$platform==platform], na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mbe_f_tmp_median <- round(sum(model_fitting_out_median$MBE_f[model_fitting_out_median$platform==platform], na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mad_d_tmp_median <- round(sum(abs(model_fitting_out_median$MBE_d[model_fitting_out_median$platform==platform]), na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mad_e_tmp_median <- round(sum(abs(model_fitting_out_median$MBE_e[model_fitting_out_median$platform==platform]), na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  mad_f_tmp_median <- round(sum(abs(model_fitting_out_median$MBE_f[model_fitting_out_median$platform==platform]), na.rm = T) / nrow(model_fitting_out_median[model_fitting_out_median$platform==platform,]),2)
  
  mbe_d_tmp_all <- round(sum(model_fitting_out_all$MBE_d[model_fitting_out_all$platform==platform], na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  mbe_e_tmp_all <- round(sum(model_fitting_out_all$MBE_e[model_fitting_out_all$platform==platform], na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  mbe_f_tmp_all <- round(sum(model_fitting_out_all$MBE_f[model_fitting_out_all$platform==platform], na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  mad_d_tmp_all <- round(sum(abs(model_fitting_out_all$MBE_d[model_fitting_out_all$platform==platform]), na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  mad_e_tmp_all <- round(sum(abs(model_fitting_out_all$MBE_e[model_fitting_out_all$platform==platform]), na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  mad_f_tmp_all <- round(sum(abs(model_fitting_out_all$MBE_f[model_fitting_out_all$platform==platform]), na.rm = T) / nrow(model_fitting_out_all[model_fitting_out_all$platform==platform,]),2)
  
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "mean",
                                       mbe_d = mbe_d_tmp_mean,
                                       mad_d = mad_d_tmp_mean,
                                       d_vals = sum(!is.na(model_fitting_out_mean$MBE_d[model_fitting_out_mean$platform==platform])),
                                       mbe_e = mbe_e_tmp_mean,
                                       mad_e = mad_e_tmp_mean,
                                       e_vals = sum(!is.na(model_fitting_out_mean$MBE_e[model_fitting_out_mean$platform==platform])),
                                       mbe_f = mbe_f_tmp_mean,
                                       mad_f = mad_f_tmp_mean,
                                       f_vals = sum(!is.na(model_fitting_out_mean$MBE_f[model_fitting_out_mean$platform==platform])),
                                       n_predictions = length(which(is.na(model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform == platform)])==F))))
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "median",
                                       mbe_d = mbe_d_tmp_median,
                                       mad_d = mad_d_tmp_median,
                                       d_vals = sum(!is.na(model_fitting_out_median$MBE_d[model_fitting_out_median$platform==platform])),
                                       mbe_e = mbe_e_tmp_median,
                                       mad_e = mad_e_tmp_median,
                                       e_vals = sum(!is.na(model_fitting_out_median$MBE_e[model_fitting_out_median$platform==platform])),
                                       mbe_f = mbe_f_tmp_median,
                                       mad_f = mad_f_tmp_median,
                                       f_vals = sum(!is.na(model_fitting_out_median$MBE_f[model_fitting_out_median$platform==platform])),
                                       n_predictions = length(which(is.na(model_fitting_out_median$SOS[which(model_fitting_out_median$platform == platform)])==F))))
  mbe_mad <- rbind(mbe_mad, data.frame(platform = platform,
                                       ndvi_mean_median = "all_values",
                                       mbe_d = mbe_d_tmp_all,
                                       mad_d = mad_d_tmp_all,
                                       d_vals = sum(!is.na(model_fitting_out_all$MBE_d[model_fitting_out_all$platform==platform])),
                                       mbe_e = mbe_e_tmp_all,
                                       mad_e = mad_e_tmp_all,
                                       e_vals = sum(!is.na(model_fitting_out_all$MBE_e[model_fitting_out_all$platform==platform])),
                                       mbe_f = mbe_f_tmp_all,
                                       mad_f = mad_f_tmp_all,
                                       f_vals = sum(!is.na(model_fitting_out_all$MBE_f[model_fitting_out_all$platform==platform])),
                                       n_predictions = length(which(is.na(model_fitting_out_all$SOS[which(model_fitting_out_all$platform == platform)])==F))))
}

mbe_d_tmp_sen1 <- round(sum(model_fitting_out_sen1$MBE_d, na.rm = T) / nrow(model_fitting_out_sen1),2)
mad_d_tmp_sen1 <- round(sum(abs(model_fitting_out_sen1$MBE_d), na.rm = T) / nrow(model_fitting_out_sen1),2)
mbe_e_tmp_sen1 <- round(sum(model_fitting_out_sen1$MBE_e, na.rm = T) / nrow(model_fitting_out_sen1),2)
mad_e_tmp_sen1 <- round(sum(abs(model_fitting_out_sen1$MBE_e), na.rm = T) / nrow(model_fitting_out_sen1),2)
mbe_f_tmp_sen1 <- round(sum(model_fitting_out_sen1$MBE_f, na.rm = T) / nrow(model_fitting_out_sen1),2)
mad_f_tmp_sen1 <- round(sum(abs(model_fitting_out_sen1$MBE_f), na.rm = T) / nrow(model_fitting_out_sen1),2)
mbe_mad <- rbind(mbe_mad, data.frame(platform = "Sentinel-1",
                                     ndvi_mean_median = NA,
                                     mbe_d = mbe_d_tmp_sen1,
                                     mad_d = mad_d_tmp_sen1,
                                     d_vals = sum(!is.na(model_fitting_out_sen1$MBE_d)),
                                     mbe_e = mbe_e_tmp_sen1,
                                     mad_e = mad_e_tmp_sen1,
                                     e_vals = sum(!is.na(model_fitting_out_sen1$MBE_e)),
                                     mbe_f = mbe_f_tmp_sen1,
                                     mad_f = mad_f_tmp_sen1,
                                     f_vals = sum(!is.na(model_fitting_out_sen1$MBE_f)),
                                     n_predictions = length(which(is.na(model_fitting_out_sen1$SOS[which(model_fitting_out_sen1$platform == "sentinel1")])==F))))


  
#perform a simple linear regression between estimated and observed budburst dates and get the rsquared of the model
lin_performance <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       ndvi_mean_median = "mean",
                                                       r_squared_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,
                                                       #slope_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[7],
                                                       r_squared_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,
                                                       #slope_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[7],
                                                       r_squared_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared))
                                                       #slope_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$coefficients[7]))
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       ndvi_mean_median = "median",
                                                       r_squared_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$r.squared,
                                                       #slope_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[7],
                                                       r_squared_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$r.squared,
                                                       #slope_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[7],
                                                       r_squared_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$r.squared))
                                                       #slope_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_median[model_fitting_out_median$platform==platform,]))$coefficients[7]))
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       ndvi_mean_median = "all_values",
                                                       r_squared_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$r.squared,
                                                       #slope_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[7],
                                                       r_squared_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$r.squared,
                                                       #slope_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[7],
                                                       r_squared_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$r.squared))
                                                       #slope_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[8],
                                                       #intercept_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_all[model_fitting_out_all$platform==platform,]))$coefficients[7]))
}
lin_performance <- rbind(lin_performance, data.frame(platform = "sentinel1",
                                                     ndvi_mean_median = NA,
                                                     r_squared_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_sen1))$r.squared,
                                                     #slope_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_sen1))$coefficients[8],
                                                     #intercept_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_sen1))$coefficients[7],
                                                     r_squared_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_sen1))$r.squared,
                                                     #slope_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_sen1))$coefficients[8],
                                                     #intercept_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_sen1))$coefficients[7],
                                                     r_squared_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_sen1))$r.squared))
                                                     #slope_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_sen1))$coefficients[8],
                                                     #intercept_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_sen1))$coefficients[7]))

lin_performance$r_squared_d <- format(round(lin_performance$r_squared_d,3), scientific=F)
lin_performance$r_squared_e <- format(round(lin_performance$r_squared_e,3), scientific=F)
lin_performance$r_squared_f <- format(round(lin_performance$r_squared_f,3), scientific=F)

# manually check for distribution of points in the fitted R² model; outliers?; few points?
#mean per-tree NDVI values
ggplot(model_fitting_out_mean, 
       aes(x=SOS, y=SOS_phase_d)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Prediction") +
  ylab("Observation") +
  ggtitle("per-tree mean NDVI values; predicted ~ Observed (Phase D), linear model")

model_fitting_out_mean %>% 
  filter(platform == "orthomosaic") %>%
  ggplot(aes(x=SOS, y=SOS_phase_d)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_light() +
  xlab("Prediction") +
  ylab("Observation") +
  ggtitle("per-tree mean NDVI values; predicted ~ Observed (Phase D), linear model")

model_fitting_out_mean %>% 
  filter(platform == "planetscope") %>%
  ggplot(aes(x=SOS, y=SOS_phase_d)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_light() +
  xlab("Prediction") +
  ylab("Observation") +
  ggtitle("per-tree mean NDVI values; predicted ~ Observed (Phase D), linear model")

model_fitting_out_mean %>% 
  filter(platform == "planetscope") %>%
  ggplot(aes(x=SOS, y=SOS_phase_d)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_light() +
  xlab("Prediction") +
  ylab("Observation") +
  ggtitle("per-tree mean NDVI values; predicted ~ Observed (Phase D), linear model")



ggplot(model_fitting_out_mean, 
       aes(x=SOS, y=SOS_phase_e)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Prediction") +
  ylab("Observation") +
  ggtitle("per-tree mean NDVI values; BB ~ observed, linear model")

ggplot(model_fitting_out_mean, 
       aes(x=SOS, y=SOS_phase_f)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Prediction") +
  ylab("Observation") +
  ggtitle("per-tree mean NDVI values; BB ~ observed, linear model")

#median per-tree NDVI values
ggplot(model_fitting_out_median, 
       aes(x=SOS, y=SOS_phase_d)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("per-tree median NDVI values; BB ~ observed, linear model")

ggplot(model_fitting_out_median, 
       aes(x=SOS, y=SOS_phase_e)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("per-tree median NDVI values; BB ~ observed, linear model")

ggplot(model_fitting_out_median, 
       aes(x=SOS, y=SOS_phase_f)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("per-tree median NDVI values; BB ~ observed, linear model")

#all per-tree NDVI values
ggplot(model_fitting_out_all, 
       aes(x=SOS, y=SOS_phase_d)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("all NDVI values; BB ~ observed, linear model")

ggplot(model_fitting_out_all, 
       aes(x=SOS, y=SOS_phase_e)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("all NDVI values; BB ~ observed, linear model")

ggplot(model_fitting_out_all, 
       aes(x=SOS, y=SOS_phase_f)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("al NDVI values; BB ~ observed, linear model")

#sentinel-1
ggplot(model_fitting_out_sen1, 
       aes(x=SOS, y=SOS_phase_d)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("Sentinel-1 backscatter; BB ~ observed, linear model")

ggplot(model_fitting_out_sen1, 
       aes(x=SOS, y=SOS_phase_e)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("Sentinel-1 backscatter; BB ~ observed, linear model")

ggplot(model_fitting_out_sen1, 
       aes(x=SOS, y=SOS_phase_f)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(platform~., scales = "free_y") +
  theme_light() +
  xlab("Day of Year") +
  ggtitle("Sentinel-1 backscatter; BB ~ observed, linear model")

#find earliest budburst, latest budburst and total number of detected budbursts per platform
platform <- "orthomosaic"
platform <- "treetalker"
platform <- "sentinel2"
platform <- "planetscope"

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

#observations
phase <- "phase_d"
phase <- "phase_e"
phase <- "phase_f"

min(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == phase)], na.rm = T)
max(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == phase)], na.rm = T)
length(which(is.na(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == phase)])==F))

