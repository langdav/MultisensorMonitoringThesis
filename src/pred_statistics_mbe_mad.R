#author: David Langenohl
#last modified: 23.08.2021
#description: derive model performance parameters; MBE (mean bias error), MAD (mean absolute deviation); R²
#NOTE:

rm(list = ls())
library(dplyr);library(ggplot2)

#load fitted models and manually observed budburst data
#load("out/log_function_models/all_values_fitted_models_output.RData")
#load("out/log_function_models/median_fitted_models_output.RData")
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")
load("out/log_function_models/budburst_fitted_models_output.RData")

model_fitting_out_mean <- rbind(model_fitting_out_mean,model_fitting_out_sen1)

#check for outliers
plot(model_fitting_out_mean$SOS[model_fitting_out_mean$platform=="planetscope"]) #no outliers

plot(model_fitting_out_mean$SOS[model_fitting_out_mean$platform=="orthomosaic"]) #far of points; none of the models was fitted right, what to do with those?!

plot(model_fitting_out_mean$SOS[model_fitting_out_mean$platform=="sentinel2"]) #remove extreme outliers; plots show, that models were not fitted accordingly
model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform=="sentinel2" & model_fitting_out_mean$tree_id %in% c("mof_cst_00007",
                                                                                                                      "mof_cst_00026",
                                                                                                                      "mof_cst_00039"))] <- NA #models not fitted right

plot(model_fitting_out_mean$SOS[model_fitting_out_mean$platform=="treetalker"]) #visually checked models that led to outliers; can be removed
model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform=="treetalker" & model_fitting_out_mean$tree_id %in% c("mof_cst_00021",
                                                                                                                       "mof_cst_00002",
                                                                                                                       "mof_cst_00012",
                                                                                                                       "mof_cst_00026"))] <- NA #models not fitted right

plot(model_fitting_out_mean$SOS[model_fitting_out_mean$platform=="sentinel1"])
model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform=="sentinel1" & model_fitting_out_mean$tree_id %in% c("mof_cst_00024",
                                                                                                                      "mof_cst_00009",
                                                                                                                      "mof_cst_00020",
                                                                                                                      "mof_cst_00001",
                                                                                                                      "mof_cst_00018",
                                                                                                                      "mof_cst_00029",
                                                                                                                      "mof_cst_00030",
                                                                                                                      "mof_cst_00040",
                                                                                                                      "mof_cst_00049"))] <- NA #models not fitted right

#testing MOS -> worse than SOS
#model_fitting_out_mean$SOS <- model_fitting_out_mean$MOS

#calculate MBE and MAD
mbe_mae <- cbind(data.frame(MBE_d = do.call(rbind, lapply(split(model_fitting_out_mean,model_fitting_out_mean$platform),function(x){round(sum(x[,"SOS"]-x[,"SOS_phase_d"],na.rm=T) / nrow(x),3)})),
                            MAE_d = do.call(rbind, lapply(split(model_fitting_out_mean,model_fitting_out_mean$platform),function(x){round(sum(abs(x[,"SOS"]-x[,"SOS_phase_d"]),na.rm=T) / nrow(x),3)})),
                            MBE_e = do.call(rbind, lapply(split(model_fitting_out_mean,model_fitting_out_mean$platform),function(x){round(sum(x[,"SOS"]-x[,"SOS_phase_e"],na.rm=T) / nrow(x),3)})),
                            MAE_e = do.call(rbind, lapply(split(model_fitting_out_mean,model_fitting_out_mean$platform),function(x){round(sum(abs(x[,"SOS"]-x[,"SOS_phase_e"]),na.rm=T) / nrow(x),3)})),
                            MBE_f = do.call(rbind, lapply(split(model_fitting_out_mean,model_fitting_out_mean$platform),function(x){round(sum(x[,"SOS"]-x[,"SOS_phase_f"],na.rm=T) / nrow(x),3)})),
                            MAE_f = do.call(rbind, lapply(split(model_fitting_out_mean,model_fitting_out_mean$platform),function(x){round(sum(abs(x[,"SOS"]-x[,"SOS_phase_f"]),na.rm=T) / nrow(x),3)}))))


#perform a simple linear regression between estimated and observed budburst dates and get the rsquared of the model
lin_performance <- NULL
for(platform in unique(model_fitting_out_mean$platform)){
  lin_performance <- rbind(lin_performance, data.frame(platform = platform,
                                                       ndvi_mean_median = "mean",
                                                       r_squared_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,
                                                       r_squared_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,
                                                       r_squared_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared))
}

lin_performance$r_squared_d <- format(round(lin_performance$r_squared_d,3), scientific=F)
lin_performance$r_squared_e <- format(round(lin_performance$r_squared_e,3), scientific=F)
lin_performance$r_squared_f <- format(round(lin_performance$r_squared_f,3), scientific=F)




mean_long %>% 
  filter(platform == "planetscope") %>%
  ggplot(aes(x=SOS, y=obs)) +
  geom_point() +
  geom_smooth(formula = y ~ x, method = "lm") +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  facet_grid(phase~., scales = "free_y") +
  theme_light() +
  xlab("Prediction") +
  ylab("Observation") +
  ggtitle("per-tree mean NDVI values; predicted ~ Observed (Phase D), linear model")

# manually check for distribution of points in the fitted R² model; outliers?; few points?
#mean per-tree NDVI values
model_fitting_out_mean <- rbind(model_fitting_out_mean, model_fitting_out_sen1)

mean_long <- rbind(cbind(model_fitting_out_mean[,1:9], obs = model_fitting_out_mean$SOS_phase_d, phase = rep("Phase D",nrow(model_fitting_out_mean))),
                   cbind(model_fitting_out_mean[,1:9], obs = model_fitting_out_mean$SOS_phase_e, phase = rep("Phase E",nrow(model_fitting_out_mean))),
                   cbind(model_fitting_out_mean[,1:9], obs = model_fitting_out_mean$SOS_phase_f, phase = rep("Phase F",nrow(model_fitting_out_mean))))

for(the_chosen_one in unique(mean_long$platform)){
  plotterich <- mean_long %>% 
    filter(platform == the_chosen_one) %>%
    filter(!is.na(SOS)) %>% 
    filter(!is.na(obs)) %>% 
    ggplot(aes(x=obs, y=SOS)) +
    geom_point() +
    geom_smooth(formula = y ~ x, method = "lm") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    facet_grid(phase~., scales = "free_y") +
    theme_light() +
    xlab("Observation DOY") +
    ylab("Prediction DOY") +
    ggtitle(paste0(the_chosen_one))
  ggsave(paste0("out/model_plots/variance_", the_chosen_one,".png"),plotterich,width = 10,height=6)
}

#find earliest budburst, latest budburst and total number of detected budbursts per platform
platform <- "orthomosaic"
platform <- "treetalker"
platform <- "sentinel2"
platform <- "planetscope"

range(model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform == platform)], na.rm = T)
length(which(is.na(model_fitting_out_mean$SOS[which(model_fitting_out_mean$platform == platform)])==F))

# range(model_fitting_out_median$SOS[which(model_fitting_out_median$platform == platform)], na.rm = T)
# length(which(is.na(model_fitting_out_median$SOS[which(model_fitting_out_median$platform == platform)])==F))
# 
# range(model_fitting_out_all$SOS[which(model_fitting_out_all$platform == platform)], na.rm = T)
# length(which(is.na(model_fitting_out_all$SOS[which(model_fitting_out_all$platform == platform)])==F))

# Sentinel-1
range(model_fitting_out_sen1$SOS, na.rm = T)
length(which(is.na(model_fitting_out_sen1$SOS)==F))

#observations
phase <- "phase_d"
phase <- "phase_e"
phase <- "phase_f"

range(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == phase)], na.rm = T)
length(which(is.na(budburst_model_fitting_out$SOS[which(budburst_model_fitting_out$phase == phase)])==F))

