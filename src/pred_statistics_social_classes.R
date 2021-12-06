#author: David Langenohl
#last modified: 06.12.2021
#description: derive model performance parameters for different social classes; MBE (mean bias error), MAD (mean absolute deviation); R²; Welch-ANOVAS
#NOTE:

rm(list = ls())
library(dplyr);library(ggplot2);library(rstatix)

#load fitted models and manually observed budburst data
load("out/log_function_models/mean_fitted_models_output.RData")
load("out/log_function_models/sentinel1_fitted_models_output.RData")

model_fitting_out_mean <- rbind(model_fitting_out_mean,model_fitting_out_sen1)
rm(model_fitting_out_sen1)

#check for (and remove) outliers
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

# load social classes data
social_classes <- read.csv("data/Data_OEP_1_Beech.csv")
social_classes <- social_classes[,c("Tree","Kraft_Class")]
colnames(social_classes) <- c("tree_id", "social_class")
social_classes <- social_classes[which(social_classes$tree_id %in% unique(model_fitting_out_mean$tree_id)),]

# loop
classes <- NULL
mbe_mae <- NULL
lin_performance <- NULL
welch_anova <- NULL
platforms <- c("orthomosaic","planetscope","sentinel1","sentinel2","treetalker")

for (class in 1:4) {
  classes <- c(classes,class)
  
  model_fitting_out_temp <- model_fitting_out_mean[which(model_fitting_out_mean$tree_id %in% unique(social_classes$tree_id[which(social_classes$social_class %in% classes)])),]
  
  #calculate MBE and MAD
  mbe_mae_tmp <- cbind(data.frame(MBE_d = do.call(rbind, lapply(split(model_fitting_out_temp,model_fitting_out_temp$platform),function(x){round(sum(x[,"SOS"]-x[,"SOS_phase_d"],na.rm=T) / nrow(x),3)})),
                                  MAE_d = do.call(rbind, lapply(split(model_fitting_out_temp,model_fitting_out_temp$platform),function(x){round(sum(abs(x[,"SOS"]-x[,"SOS_phase_d"]),na.rm=T) / nrow(x),3)})),
                                  MBE_e = do.call(rbind, lapply(split(model_fitting_out_temp,model_fitting_out_temp$platform),function(x){round(sum(x[,"SOS"]-x[,"SOS_phase_e"],na.rm=T) / nrow(x),3)})),
                                  MAE_e = do.call(rbind, lapply(split(model_fitting_out_temp,model_fitting_out_temp$platform),function(x){round(sum(abs(x[,"SOS"]-x[,"SOS_phase_e"]),na.rm=T) / nrow(x),3)})),
                                  MBE_f = do.call(rbind, lapply(split(model_fitting_out_temp,model_fitting_out_temp$platform),function(x){round(sum(x[,"SOS"]-x[,"SOS_phase_f"],na.rm=T) / nrow(x),3)})),
                                  MAE_f = do.call(rbind, lapply(split(model_fitting_out_temp,model_fitting_out_temp$platform),function(x){round(sum(abs(x[,"SOS"]-x[,"SOS_phase_f"]),na.rm=T) / nrow(x),3)}))))
  
  mbe_mae_tmp$social_classes <- paste0(classes, sep = "", collapse = "")
  mbe_mae_tmp$platform <- platforms
  mbe_mae_tmp <- mbe_mae_tmp[,c(7,8,1:6)]
  mbe_mae <- rbind(mbe_mae, mbe_mae_tmp)
  rownames(mbe_mae_tmp) <- NULL
  
  #perform a simple linear regression between estimated and observed budburst dates and get the rsquared of the model
  lin_performance_tmp <- NULL
  for(platform in platforms){
    lin_performance_tmp <- rbind(lin_performance_tmp, data.frame(platform = platform,
                                                                 ndvi_mean_median = "mean",
                                                                 r_squared_d = format(round(summary(lm(SOS ~ SOS_phase_d,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,3), scientific=F),
                                                                 r_squared_e = format(round(summary(lm(SOS ~ SOS_phase_e,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,3), scientific=F),
                                                                 r_squared_f = format(round(summary(lm(SOS ~ SOS_phase_f,model_fitting_out_mean[model_fitting_out_mean$platform==platform,]))$r.squared,3), scientific=F)))
  }
  
  lin_performance_tmp$social_classes <- paste0(classes, sep = "", collapse = "")
  lin_performance_tmp <- lin_performance_tmp[,c(6,1:5)]
  lin_performance <- rbind(lin_performance, lin_performance_tmp)
  
  #Welch-Anova(s) to test, if per-platform estimated values differ from observed values
  welch_anova_tmp <- NULL
  for(platform in platforms){
    
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
    
    
    welch_anova_tmp <- rbind(welch_anova_tmp, data.frame(platform = platform,
                                                         mean_median = "mean_d",
                                                         p_val = rstatix::welch_anova_test(tmp_df_mean_d, values ~ estimated_observed)$p,
                                                         f = as.numeric(rstatix::welch_anova_test(tmp_df_mean_d, values ~ estimated_observed)$statistic)))
    
    
    welch_anova_tmp <- rbind(welch_anova_tmp, data.frame(platform = platform,
                                                         mean_median = "mean_e",
                                                         p_val = rstatix::welch_anova_test(tmp_df_mean_e, values ~ estimated_observed)$p,
                                                         f = as.numeric(rstatix::welch_anova_test(tmp_df_mean_e, values ~ estimated_observed)$statistic)))
    
    welch_anova_tmp <- rbind(welch_anova_tmp, data.frame(platform = platform,
                                                         mean_median = "mean_f",
                                                         p_val = rstatix::welch_anova_test(tmp_df_mean_f, values ~ estimated_observed)$p,
                                                         f = as.numeric(rstatix::welch_anova_test(tmp_df_mean_f, values ~ estimated_observed)$statistic)))
    
  }
  
  welch_anova_tmp$p_val <- format(round(welch_anova_tmp$p_val,6), scientific=F)
  welch_anova_tmp$social_classes <- paste0(classes, sep = "", collapse = "")
  welch_anova_tmp <- welch_anova_tmp[,c(5,1:4)]
  welch_anova <- rbind(welch_anova, welch_anova_tmp)
  
}
rm(lin_performance_tmp);rm(mbe_mae_tmp);rm(model_fitting_out_temp);rm(tmp_df_mean_d);rm(tmp_df_mean_e);rm(tmp_df_mean_f);rm(welch_anova_tmp)
