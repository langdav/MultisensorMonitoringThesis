#author: David Langenohl; based on the work of Dominic Fawcett
#last modified: 19.08.2021
#description: script to extract NDVI time series per crown and fit logistic function to retrieve phenological start and end of budburst phase
#NOTE: no models fitted to Sentinel-2 data, as too little data available (3 days)
#NOTE: model variables: a = amplitude of increase; b = growth rate; c = value of t at inflection point; d = lower asymptote

rm(list = ls())
library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);library(kit);library(minpack.lm);library(rstatix) #minpack.lm isneeded for the use of minpack.lm::nlsLM()

#load data
load("out/all_in_one/aio_daily_ndvi_all.RData")
aio_all <- aio_daily_ndvi_all;rm(aio_daily_ndvi_all)

multiplatform <- aio_all %>% 
  filter(platform %in% c("planetscope","sentinel2","orthomosaic")) %>% 
  group_by(tree_id, date) %>% 
  dplyr::summarize(ndvi_mean = mean(ndvi, na.rm = T),
                   ndvi_sd = sd(ndvi, na.rm = T),
                   budburst_perc = mean(budburst_perc, na.rm = T)) %>% 
  ungroup()

load("data/trees.RData")
trees$tree_id <- as.character(trees$tree_id)
trees <- trees[c(1:50),] #reduce to trees with a budburst record

#model fitting 
model_fitting_out <- NULL
models <- list()
countvar <- 1

#starting values
a_start <- 0.5    #lower asymptote
b_start <- 1      #amplitude of NDVI
c_start <- 0.7    #inflection point
d_start <- 0.01   #growth rate

for(tree in unique(multiplatform$tree_id)){
  
  #get doys; as not all trees are present in some of the orthomosaics, the doylist is extracted tree-specific
  doylist <- sort(yday(multiplatform$date[multiplatform$tree_id == tree]))
  doystart <- head(doylist,1)
  doyend <- tail(doylist,1)
  
  ndvidat <- multiplatform$ndvi_mean[which(multiplatform$tree_id == tree)]
  sd <- multiplatform$ndvi_sd[which(multiplatform$tree_id == tree)]
  
  if(length(ndvidat) == 0){
    print(paste('not data for',', tree ',tree))
    model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "multiplatform",
                                                             tree_id = tree,
                                                             SOS = NA,
                                                             MOS = NA,
                                                             EOS = NA,
                                                             RSE = NA,
                                                             no_data = T,
                                                             warning = F,
                                                             error = F)) #residualSTDerror
    
    models[[countvar]] <- list(platform = "multiplatform",
                               tree_id = tree,
                               ndvidat = NA,
                               ndvi_sd = sd,
                               doylist = NA,
                               model = NA)
    countvar <- countvar+1
  } else {
    
    # see, if model can be calculated; if not, move on
    # tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
    tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ b/(1 + exp(c + (doylist*d))) + a,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
    
    if(is(tt,'warning')){
      print(paste('warning for', 'tree ',tree))
      model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "multiplatform",
                                                               tree_id = tree,
                                                               SOS = NA,
                                                               MOS = NA,
                                                               EOS = NA,
                                                               RSE = NA,
                                                               no_data = F,
                                                               warning = T,
                                                               error = F)) #residualSTDerror
      
      models[[countvar]] <- list(platform = "multiplatform", 
                                 tree_id = tree,
                                 ndvidat = ndvidat,
                                 ndvi_sd = sd,
                                 doylist = doylist, 
                                 model = NA)
      countvar <- countvar+1
    } else if(is(tt,'error')){
      print(paste0('error at platform ',"multiplatform", ', tree ',tree,', ',tt))
      model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "multiplatform",
                                                               tree_id = tree,
                                                               SOS = NA,
                                                               MOS = NA,
                                                               EOS = NA,
                                                               RSE = NA,
                                                               no_data = F,
                                                               warning = F,
                                                               error = T)) #residualSTDerror
      
      models[[countvar]] <- list(platform = "multiplatform", 
                                 tree_id = tree, 
                                 ndvidat = ndvidat,
                                 ndvi_sd = sd,
                                 doylist = doylist, 
                                 model = NA)
      countvar <- countvar+1
    } else {
      # minpack.lm isneeded for the use of minpack.lm::nlsLM() instead of stats::nls(), as it is more robust to bad starting values (and I couldn't find good ones)
      # it uses the Levenberg-Marquardt (https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) algorithm instead of Gauss-Newton
      fitmodel <- minpack.lm::nlsLM(ndvidat ~ b/(1 + exp(c + (doylist*d))) + a,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start))
      
      #plot(ndvidat ~ doylist)
      #curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)
      
      #get model coefficients
      modsum <- summary(fitmodel) #all parameters significant
      a <- modsum$parameters[1]
      b <- modsum$parameters[2]
      c <- modsum$parameters[3]
      d <- modsum$parameters[4]
      
      # predict values of days, that are not present in the images, based on lgistic function "fitmodel"
      preddoylist <- seq(doystart,doyend,1)
      predNDVImod <- predict(fitmodel,data.frame(doylist=preddoylist))
      
      #curvature change rate (= first derivation)
      firstDerivCurv <- -1*(b*(d^3)*exp(c+(preddoylist*d))*(exp(6*c+(6*preddoylist*d))+(-2*(b^2)*exp(4*c)*(d^2)-9*exp(4*c))*exp(4*d*preddoylist)+(2*(b^2)*exp(3*c)*(d^2)-16*exp(3*c))*exp(3*d*preddoylist)+(-2*(b^2)*exp(2*c)*(d^2)-9*exp(2*c))*exp(2*d*preddoylist)+1) / (((exp(c+(preddoylist*d))+1)^8)*((((((b^2)*(d^2)*exp(2*(c+(preddoylist*d)))) / ((exp(c+(preddoylist*d))+1)^4))) + 1)^2.5))) #new model
      
      #plot
      # ggplot() +
      #   geom_point(aes(preddoylist,predNDVImod), colour = "red") +
      #   geom_point(aes(doylist,ndvidat)) +
      #   scale_x_continuous(breaks = round(seq(min(doylist), max(doylist), by = 1),1)) +
      #   scale_y_continuous(breaks = seq(0, 100, by = 5)) +
      #   geom_line(aes(preddoylist,firstDerivNDVImod),colour = "red") +
      #   geom_line(aes(preddoylist,secondDerivNDVImod),colour = "blue") +
      #   geom_line(aes(preddoylist,firstDerivCurv*10),colour = "red") +
      #   geom_line(aes(preddoylist, ROCcurvature*10),colour = "green") +
      #   geom_vline(xintercept = preddoylist[which(secondDerivNDVImod == max(secondDerivNDVImod))]) +
      #   geom_vline(xintercept = preddoylist[which(secondDerivNDVImod == min(secondDerivNDVImod))]) +
      #   geom_hline(yintercept = predNDVImod[which(secondDerivNDVImod == max(secondDerivNDVImod))]) +
      #   geom_hline(yintercept = predNDVImod[which(secondDerivNDVImod == min(secondDerivNDVImod))])
      
      
      if(any(is.na(firstDerivCurv)==T)){
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "multiplatform",
                                                                 tree_id = tree,
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = modsum$sigma,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
      } else {
        # new way of finding SOS
        # SOS = first max. of firstDerivCurv; MOS = min. of firstDerivCurv; EOS = second max. of firstDerivCurv; 
        SOS <- sort(preddoylist[kit::topn(firstDerivCurv,2)])[1]
        MOS <- preddoylist[which.min(firstDerivCurv)]
        EOS <- sort(preddoylist[kit::topn(firstDerivCurv,2)])[2]
        
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "multiplatform",
                                                                 tree_id = tree,
                                                                 SOS = SOS,
                                                                 MOS = MOS,
                                                                 EOS = EOS,
                                                                 RSE = modsum$sigma,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
        # example plot for thesis
        # ggplot() +
        #   geom_line(aes(preddoylist,predNDVImod)) +
        #   geom_point(aes(doylist,ndvidat)) +
        #   geom_errorbar(aes(x = doylist ,ymin=ndvidat-sd, ymax=ndvidat+sd), width = 1) +
        #   geom_line(aes(preddoylist,firstDerivCurv*1000+0.7), linetype = "dashed") +
        #   geom_vline(xintercept = SOS, linetype = "longdash", colour = "red") +
        #   geom_text(aes(x = SOS-1,
        #                 y = 0.4,
        #                 label = "budburst",
        #                 vjust = 0,
        #                 angle = 90)) +
        #   theme_minimal() +
        #   xlab("DOY") +
        #   ylab("NDVI")
        # ggsave(paste0("tmp/example_curvature_derivative.png"),width = 5,height=4)
      }
      
      models[[countvar]] <- list(platform = "multiplatform", 
                                 tree_id = tree,
                                 ndvidat = ndvidat,
                                 ndvi_sd = sd,
                                 doylist = doylist,
                                 model = fitmodel)
      
      countvar <- countvar+1
    }
  }
}

model_fitting_out_multiplatform <- model_fitting_out
model_out_multiplatform <- models


#add predicted budburst from observations
load("out/log_function_models/budburst_fitted_models_output.RData")
obs_budburst <- merge(budburst_model_fitting_out[which(budburst_model_fitting_out$phase == "phase_d"),c(2,3)],
                      budburst_model_fitting_out[which(budburst_model_fitting_out$phase == "phase_e"),c(2,3)],by="tree_id")
obs_budburst <- merge(obs_budburst, budburst_model_fitting_out[which(budburst_model_fitting_out$phase == "phase_f"),c(2,3)],by="tree_id")
colnames(obs_budburst) <- c("tree_id", "SOS_phase_d", "SOS_phase_e", "SOS_phase_f")

model_fitting_out_multiplatform <- merge(model_fitting_out_multiplatform, obs_budburst, by = "tree_id", all.x =T)

save(model_fitting_out_multiplatform, file = "out/log_function_models/multiplatform_fitted_models_output.RData")
save(model_out_multiplatform, file = "out/log_function_models/multiplatforms_fitted_models.RData")

#look at some of those models
i <- 1
ndvidat <- model_out_multiplatform[[i]]$ndvidat;doylist <- model_out_multiplatform[[i]]$doylist;fitmodel <- model_out_multiplatform[[i]]$model;plot(ndvidat ~ doylist, main = model_out_multiplatform[[i]]$platform);curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)  

#plotting
plot_multi <- function(i){
  
  ndvidat <- model_out_multiplatform[[i]]$ndvidat
  doylist <- model_out_multiplatform[[i]]$doylist
  fitmodel <- model_out_multiplatform[[i]]$model
  ndvi_sd <- model_out_multiplatform[[i]]$ndvi_sd
  platform <- "multiplatform"
  tree <- model_out_multiplatform[[i]]$tree
  preddoylist <- seq(head(doylist, 1),tail(doylist,1),1)
  predNDVImod <- predict(fitmodel,data.frame(doylist=preddoylist))
  
  plot <- ggplot() +
    geom_point(aes(doylist, ndvidat)) +
    geom_errorbar(aes(x = doylist ,ymin=ndvidat-ndvi_sd, ymax=ndvidat+ndvi_sd), width = 1) +
    geom_line(aes(preddoylist,predNDVImod)) +
    geom_vline(xintercept = model_fitting_out_multiplatform$SOS[i], linetype = "dotdash", color = "blue", size = .8) +
    geom_text(aes(x = model_fitting_out_multiplatform$SOS[i],
                  y = min(ndvidat),
                  label = "predicted",
                  vjust = 0,
                  angle = 90)) +
    geom_vline(xintercept = model_fitting_out_multiplatform$SOS_phase_d[i], linetype = "longdash", color = "red", size = .5) +
    geom_text(aes(x = model_fitting_out_multiplatform$SOS_phase_d[i],
                  y = max(ndvidat),
                  label = "phase_d",
                  vjust = 0,
                  angle = 90)) +
    geom_vline(xintercept = model_fitting_out_multiplatform$SOS_phase_e[i], linetype = "longdash", color = "red", size = .5) +
    geom_text(aes(x = model_fitting_out_multiplatform$SOS_phase_e[i],
                  y = max(ndvidat),
                  label = "phase_e",
                  vjust = 0,
                  angle = 90)) +
    geom_vline(xintercept = model_fitting_out_multiplatform$SOS_phase_f[i], linetype = "longdash", color = "red", size = .5) +
    geom_text(aes(x = model_fitting_out_multiplatform$SOS_phase_f[i],
                  y = max(ndvidat),
                  label = "phase_f",
                  vjust = 0,
                  angle = 90)) +
    theme_light() +
    xlab("Day of Year") +
    ylab("NDVI") +
    ggtitle(paste0("Multiplatform; ",tree))
  return(plot)
}

for(i in seq(1,37,12)){
  v <- i+11
  
  all_trees <- lapply(i:v, function(x) plot_multi(x))
  all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 4, ncol = 3,layout_matrix = matrix(1:12, 4, 3, TRUE))
  ggsave(paste0("out/model_plots/multiplatform_",i,"_",v,".png"),all_trees_panel,width = 16,height=16)
}

all_trees <- lapply(49:50, function(x) plot_multi(x))
all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 1, ncol = 2,layout_matrix = matrix(1:2, 1, 2, TRUE))
ggsave(paste0("out/model_plots/multiplatform_",48,"_",50,".png"),all_trees_panel,width = 10.7,height=4)

#remove outliers
plot(model_fitting_out_multiplatform$SOS) #one potential outlier; after checking it visually, it can't be removed


#ealierst/latest budburst dates
range(model_fitting_out_multiplatform$SOS, na.rm = T)
length(which(is.na(model_fitting_out_multiplatform$SOS)==F))

#statistics
#calculate MBE and MAD
model_fitting_out_multiplatform$MBE_d <- model_fitting_out_multiplatform$SOS - model_fitting_out_multiplatform$SOS_phase_d
model_fitting_out_multiplatform$MBE_e <- model_fitting_out_multiplatform$SOS - model_fitting_out_multiplatform$SOS_phase_e
model_fitting_out_multiplatform$MBE_f <- model_fitting_out_multiplatform$SOS - model_fitting_out_multiplatform$SOS_phase_f

mbe_mad <- data.frame(platform = "multiplatform",
                      mbe_d = round(sum(model_fitting_out_multiplatform$MBE_d, na.rm = T) / nrow(model_fitting_out_multiplatform),2),
                      mad_d = round(sum(abs(model_fitting_out_multiplatform$MBE_d), na.rm = T) / nrow(model_fitting_out_multiplatform),2),
                      d_vals = sum(!is.na(model_fitting_out_multiplatform$MBE_d)),
                      mbe_e = round(sum(model_fitting_out_multiplatform$MBE_e, na.rm = T) / nrow(model_fitting_out_multiplatform),2),
                      mad_e = round(sum(abs(model_fitting_out_multiplatform$MBE_e), na.rm = T) / nrow(model_fitting_out_multiplatform),2),
                      e_vals = sum(!is.na(model_fitting_out_multiplatform$MBE_e)),
                      mbe_f = round(sum(model_fitting_out_multiplatform$MBE_f, na.rm = T) / nrow(model_fitting_out_multiplatform),2),
                      mad_f = round(sum(abs(model_fitting_out_multiplatform$MBE_f), na.rm = T) / nrow(model_fitting_out_multiplatform),2),
                      f_vals = sum(!is.na(model_fitting_out_multiplatform$MBE_f)),
                      n_predictions = length(which(is.na(model_fitting_out_multiplatform$SOS)==F)))

#perform a simple linear regression between estimated and observed budburst dates and get the rsquared of the model
lin_performance <- data.frame(platform = "multiplatform",
                              r_squared_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_multiplatform))$r.squared,
                              # slope_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_multiplatform))$coefficients[8],
                              # intercept_signif_d = summary(lm(SOS ~ SOS_phase_d,model_fitting_out_multiplatform))$coefficients[7],
                              r_squared_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_multiplatform))$r.squared,
                              # slope_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_multiplatform))$coefficients[8],
                              # intercept_signif_e = summary(lm(SOS ~ SOS_phase_e,model_fitting_out_multiplatform))$coefficients[7],
                              r_squared_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_multiplatform))$r.squared)
# slope_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_multiplatform))$coefficients[8],
# intercept_signif_f = summary(lm(SOS ~ SOS_phase_f,model_fitting_out_multiplatform))$coefficients[7])

multi_long <- rbind(cbind(model_fitting_out_multiplatform[,1:9], obs = model_fitting_out_multiplatform$SOS_phase_d, phase = rep("Phase D",nrow(model_fitting_out_multiplatform))),
                    cbind(model_fitting_out_multiplatform[,1:9], obs = model_fitting_out_multiplatform$SOS_phase_e, phase = rep("Phase E",nrow(model_fitting_out_multiplatform))),
                    cbind(model_fitting_out_multiplatform[,1:9], obs = model_fitting_out_multiplatform$SOS_phase_f, phase = rep("Phase F",nrow(model_fitting_out_multiplatform))))


plotterich <- multi_long %>% 
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
  ggtitle("multiplatform")
ggsave(paste0("out/model_plots/variance_multiplatform.png"),plotterich,width = 10,height=6)


#ANOVAS
tmp_df_all_d <- data.frame(values = c(model_fitting_out_multiplatform$SOS,
                                      model_fitting_out_multiplatform$SOS_phase_d),
                           estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_multiplatform)),
                                                  rep(as.factor("observed"), nrow(model_fitting_out_multiplatform))))

tmp_df_all_e <- data.frame(values = c(model_fitting_out_multiplatform$SOS,
                                      model_fitting_out_multiplatform$SOS_phase_e),
                           estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_multiplatform)),
                                                  rep(as.factor("observed"), nrow(model_fitting_out_multiplatform))))

tmp_df_all_f <- data.frame(values = c(model_fitting_out_multiplatform$SOS,
                                      model_fitting_out_multiplatform$SOS_phase_f),
                           estimated_observed = c(rep(as.factor("estimated"), nrow(model_fitting_out_multiplatform)),
                                                  rep(as.factor("observed"), nrow(model_fitting_out_multiplatform))))

rstatix::welch_anova_test(tmp_df_all_d, values ~ estimated_observed)$p
as.numeric(rstatix::welch_anova_test(tmp_df_all_d, values ~ estimated_observed)$statistic)
rstatix::welch_anova_test(tmp_df_all_e, values ~ estimated_observed)$p
as.numeric(rstatix::welch_anova_test(tmp_df_all_e, values ~ estimated_observed)$statistic)
rstatix::welch_anova_test(tmp_df_all_f, values ~ estimated_observed)$p
as.numeric(rstatix::welch_anova_test(tmp_df_all_f, values ~ estimated_observed)$statistic)
