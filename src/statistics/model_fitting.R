#author: David Langenohl; based on the work of Dominic Fawcett
#last modified: 19.08.2021
#description: script to extract NDVI time series per crown and fit logistic function to retrieve phenological start and end of budburst phase
#NOTE: no models fitted to Sentinel-2 data, as too little data available (3 days)
#NOTE: model variables: a = amplitude of increase; b = growth rate; c = value of t at inflection point; d = lower asymptote

rm(list = ls())
library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);library(minpack.lm) #minpack.lm isneeded for the use of minpack.lm::nlsLM()

#load data
load("out/all_in_one/aio_daily_ndvi_all.RData")
aio_all <- aio_daily_ndvi_all;rm(aio_daily_ndvi_all)
aio_all <- aio_all[-which(aio_all$date > as.Date("2021-07-01")),] #limit to data before "2021-07-01"

load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
aio_mean <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
aio_mean <- aio_mean[-which(aio_mean$date > as.Date("2021-07-01")),] #limit to data before "2021-07-01"

load("out/all_in_one/aio_daily_ndvi_per_tree_medians.RData")
aio_median <- aio_daily_ndvi_medians;rm(aio_daily_ndvi_medians)
aio_median <- aio_median[-which(aio_median$date > as.Date("2021-07-01")),] #limit to data before "2021-07-01"

load("out/sentinel1/backscatter_all_with_phenoclasses_sentinel1.RData")

#model fitting function
model_fitting <- function(dataset = aio_mean, mean = T, median = F, return_models = F, sentinel1 = F){
  #perform model fitting and extract SOS, MOS and EOS values from fitted model
  model_fitting_out <- NULL
  models <- list()
  countvar <- 1
  
  aio <- dataset
  
  if(sentinel1 == F){
    #iterate through platforms
    for(platform in unique(aio$platform)){
      
      SOSdoy <- list()
      platform_data_only <- aio[aio$platform == platform,]
      
      ntrees = 50
      SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)
      
      #starting values
      # a_start <- 0.6    #date of budburst
      # b_start <- 0.01   #growth rate
      # c_start <- 1    #amplitude of increase
      # d_start <- 0    #lower asymptote
      
      # if(platform == "planetscope"){
      #   a_start <- 0.6    #date of budburst
      #   b_start <- 0.01   #growth rate
      #   c_start <- 0.5    #amplitude of increase
      #   d_start <- 0.5    #lower asymptote
      # } else if(platform == "sentinel2"){
      #   a_start <- 120
      #   b_start <- 0.01
      #   c_start <- 0.3
      #   d_start <- 0.5
      # } else if(platform == "treetalker"){
      #   a_start <- 120
      #   b_start <- 0.01
      #   c_start <- 0.8
      #   d_start <- 0
      # } else if(platform == "orthomosaic"){
      #   a_start <- 120
      #   b_start <- 0.001
      #   c_start <- 0.8
      #   d_start <- 0
      # }
      
      
      if(platform == "planetscope"){
        a_start <- 0.5    #lower asymptote
        b_start <- 1      #amplitude of NDVI
        c_start <- 0.7    #inflection point
        d_start <- 0.01   #growth rate
      } else if(platform == "sentinel2"){
        a_start <- 0.5    #lower asymptote
        b_start <- 1      #amplitude of NDVI
        c_start <- 0.7    #inflection point
        d_start <- 0.01   #growth rate
      } else if(platform == "treetalker"){
        a_start <- 0      #lower asymptote
        b_start <- 1      #amplitude of NDVI
        c_start <- 0.7    #inflection point
        d_start <- 0.01   #growth rate
      } else if(platform == "orthomosaic"){
        a_start <- 0      #lower asymptote
        b_start <- 1      #amplitude of NDVI
        c_start <- 0.7    #inflection point
        d_start <- 0.01   #growth rate
      }
      
      #iterate through all trees per platform
      for(tree in unique(aio$tree_id)){
        
        #get doys; as not all trees are present in some of the orthomosaics, the doylist is extracted tree-specific
        #doylist <- sort(yday(unique(platform_data_only$date[platform_data_only$tree_id == tree])))
        doylist <- sort(yday(platform_data_only$date[platform_data_only$tree_id == tree]))
        doystart <- head(doylist,1)
        doyend <- tail(doylist,1)
        
        if(mean == T & median == F){
          ndvidat <- platform_data_only$ndvi_mean[which(platform_data_only$tree_id == tree)]
        } else if(mean == F & median == T){
          ndvidat <- platform_data_only$ndvi_median[which(platform_data_only$tree_id == tree)]
        } else {
          ndvidat <- platform_data_only$ndvi[which(platform_data_only$tree_id == tree)]
        }
        # ifelse(mean == T,
        #        ndvidat <- platform_data_only$ndvi_mean[which(platform_data_only$tree_id == tree)],
        #        ndvidat <- platform_data_only$ndvi_median[which(platform_data_only$tree_id == tree)])
        
        
        if(length(ndvidat) == 0){
          print(paste('not data for platform ',platform, ', tree ',tree))
          model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                   tree_id = tree,
                                                                   SOS = NA,
                                                                   MOS = NA,
                                                                   EOS = NA,
                                                                   RSE = NA,
                                                                   no_data = T,
                                                                   warning = F,
                                                                   error = F)) #residualSTDerror
          
          models[[countvar]] <- list(platform = platform, tree_id = tree, ndvidat = NA, doylist = NA, model = NA)
          countvar <- countvar+1
        } else {
          
          # see, if model can be calculated; if not, move on
          # tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
          tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ b/(1 + exp(c + (doylist*d))) + a,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
          
          if(is(tt,'warning')){
            print(paste('warning at platform ',platform, ', tree ',tree))
            model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                     tree_id = tree,
                                                                     SOS = NA,
                                                                     MOS = NA,
                                                                     EOS = NA,
                                                                     RSE = NA,
                                                                     no_data = F,
                                                                     warning = T,
                                                                     error = F)) #residualSTDerror
            
            models[[countvar]] <- list(platform = platform, tree_id = tree, ndvidat = ndvidat, doylist = doylist, model = NA)
            countvar <- countvar+1
          } else if(is(tt,'error')){
            print(paste0('error at platform ',platform, ', tree ',tree,', ',tt))
            model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                     tree_id = tree,
                                                                     SOS = NA,
                                                                     MOS = NA,
                                                                     EOS = NA,
                                                                     RSE = NA,
                                                                     no_data = F,
                                                                     warning = F,
                                                                     error = T)) #residualSTDerror
            
            models[[countvar]] <- list(platform = platform, tree_id = tree, ndvidat = ndvidat, doylist = doylist, model = NA)
            countvar <- countvar+1
          } else {
            # minpack.lm isneeded for the use of minpack.lm::nlsLM() instead of stats::nls(), as it is more robust to bad starting values (and I couldn't find good ones)
            # it uses the Levenberg-Marquardt (https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) algorithm instead of Gauss-Newton
            # fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start))
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
            predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
            
            #first derivation: slope
            firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2
            
            #second derivation: curvature
            secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2))
            curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
            
            #change of curvature
            ROCcurvature <- c(curvature[2:(doyend-doystart+1)],NA)-curvature 
            
            # diff(ROCcurvature) = amount of change of curvature between two points
            # sign(diff(ROCcurvature)) = positive (1) or negative (-1) change of curvature between two points
            # which(diff(sign(diff(ROCcurvature)))==-2) = point where curvature changes from positive to negative (or the other way around)
            SOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#get DOY of ROC local maximum 1
            MOS <- preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]
            EOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[2]]#get DOY of ROC local maximum 2
            
            model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                     tree_id = tree,
                                                                     SOS = SOS,
                                                                     MOS = MOS,
                                                                     EOS = EOS,
                                                                     RSE = modsum$sigma,
                                                                     no_data = F,
                                                                     warning = F,
                                                                     error = F)) #residualSTDerror
            
            models[[countvar]] <- list(platform = platform, tree_id = tree, ndvidat = ndvidat, doylist = doylist, model = fitmodel)
            countvar <- countvar+1
          }
        }
      }
    }
  } else {
    
    SOSdoy <- list()
    ntrees = 50
    SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)
    
    #starting values
    a_start <- 5    #lower asymptote
    b_start <- 10      #amplitude of backscatter
    c_start <- 6    #inflection point
    d_start <- 0.01   #growth rate
    
    #iterate through all trees per platform
    for(tree in unique(aio$tree_id)){
      
      #get doys; as not all trees are present in some of the orthomosaics, the doylist is extracted tree-specific
      doylist <- sort(aio$doy[aio$tree_id == tree])
      doystart <- head(doylist,1)
      doyend <- tail(doylist,1)
      
      ndvidat <- aio$sigma_ratio[which(aio$tree_id == tree)]
      
      if(length(ndvidat) == 0){
        print(paste('not data for platform ',platform, ', tree ',tree))
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                 tree_id = tree,
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = NA,
                                                                 no_data = T,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
        
        models[[countvar]] <- list(platform = "sentinel1", tree_id = tree, ndvidat = NA, doylist = NA, model = NA)
        countvar <- countvar+1
      } else {
        # see, if model can be calculated; if not, move on
        # tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
        tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ b/(1 + exp(c + (doylist*d))) + a,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
        
        if(is(tt,'warning')){
          print(paste('warning at platform ',"sentinel1", ', tree ',tree))
          model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                   tree_id = tree,
                                                                   SOS = NA,
                                                                   MOS = NA,
                                                                   EOS = NA,
                                                                   RSE = NA,
                                                                   no_data = F,
                                                                   warning = T,
                                                                   error = F)) #residualSTDerror
          
          models[[countvar]] <- list(platform = "sentinel1", tree_id = tree, ndvidat = ndvidat, doylist = doylist, model = NA)
          countvar <- countvar+1
        } else if(is(tt,'error')){
          print(paste('error at platform ',"sentinel1", ', tree ',tree))
          model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                   tree_id = tree,
                                                                   SOS = NA,
                                                                   MOS = NA,
                                                                   EOS = NA,
                                                                   RSE = NA,
                                                                   no_data = F,
                                                                   warning = F,
                                                                   error = T)) #residualSTDerror
          
          models[[countvar]] <- list(platform = "sentinel1", tree_id = tree, ndvidat = ndvidat, doylist = doylist, model = NA)
          countvar <- countvar+1
        } else {
          # minpack.lm isneeded for the use of minpack.lm::nlsLM() instead of stats::nls(), as it is more robust to bad starting values (and I couldn't find good ones)
          # it uses the Levenberg-Marquardt (https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) algorithm instead of Gauss-Newton
          # fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start))
          fitmodel <- minpack.lm::nlsLM(ndvidat ~ b/(1 + exp(c + (doylist*d))) + a,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start))
          
          #plot(ndvidat ~ doylist)
          #curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)
          
          #get model coefficients
          modsum <- summary(fitmodel)
          a <- modsum$parameters[1]
          b <- modsum$parameters[2]
          c <- modsum$parameters[3]
          d <- modsum$parameters[4]
          
          # predict values of days, that are not present in the images, based on lgistic function "fitmodel"
          preddoylist <- seq(doystart,doyend,1)
          predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
          
          #first derivation: slope
          firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2
          
          #second derivation: curvature
          secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2))
          curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
          
          #change of curvature
          ROCcurvature <- c(curvature[2:(doyend-doystart+1)],NA)-curvature 
          
          # diff(ROCcurvature) = amount of change of curvature between two points
          # sign(diff(ROCcurvature)) = positive (1) or negative (-1) change of curvature between two points
          # which(diff(sign(diff(ROCcurvature)))==-2) = point where curvature changes from positive to negative (or the other way around)
          SOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#get DOY of ROC local maximum 1
          MOS <- preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]
          EOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[2]]#get DOY of ROC local maximum 2
          
          model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                   tree_id = tree,
                                                                   SOS = SOS,
                                                                   MOS = MOS,
                                                                   EOS = EOS,
                                                                   RSE = modsum$sigma,
                                                                   no_data = F,
                                                                   warning = F,
                                                                   error = F)) #residualSTDerror
          
          models[[countvar]] <- list(platform = "sentinel1", tree_id = tree, ndvidat = ndvidat, doylist = doylist, model = fitmodel)
          countvar <- countvar+1
        }
      }
    } 
  }
  model_fitting_out <- model_fitting_out[!duplicated(model_fitting_out[,c(1,2)]),]
  ifelse(return_models == T,return(models),return(model_fitting_out))
}

model_fitting_out_all <- model_fitting(dataset = aio_all, mean = F, median = F, return_models = F, sentinel1 = F)
model_fitting_out_mean <- model_fitting(dataset = aio_mean, mean = T, median = F, return_models = F, sentinel1 = F)
model_fitting_out_median <- model_fitting(dataset = aio_median, mean = F, median = T,return_models = F, sentinel1 = F)

models_all <- model_fitting(dataset = aio_all, mean = F, median = F, return_models = T, sentinel1 = F)
models_mean <- model_fitting(dataset = aio_mean, mean = T, median = F, return_models = T, sentinel1 = F)
models_median <- model_fitting(dataset = aio_median, mean = F, median = T, return_models = T, sentinel1 = F)

model_fitting_out_sen1 <- model_fitting(dataset = sen1_backscatter, mean = T, return_models = F, sentinel1 = T)
models_sen1 <- model_fitting(dataset = sen1_backscatter, mean = T, median = F,return_models = T, sentinel1 = T)

#look at some of those models
i <- 45
ndvidat <- models_mean[[i]]$ndvidat;doylist <- models_mean[[i]]$doylist;fitmodel <- models_mean[[i]]$model;plot(ndvidat ~ doylist, main = models_mean[[i]]$platform);curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)  

#sentinel-1
i <- 2
ndvidat <- models_sen1[[i]]$ndvidat;doylist <- models_sen1[[i]]$doylist;fitmodel <- models_sen1[[i]]$model;plot(ndvidat ~ doylist, main = models_sen1[[i]]$platform);curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)  



#for each tree, add doy of manually observed budburst
budburst <- read.csv("data/budburst_data/budburst_long.csv")
budburst$budburst_obervation_doy <- yday(budburst$date)
budburst <- budburst[which(budburst$budburst==T),c("tree_id","budburst_obervation_doy","budburst_perc")]
budburst <- budburst[!duplicated(budburst$tree_id),]

model_fitting_out_all <- merge(model_fitting_out_all, budburst, by = "tree_id", all.x =T)
model_fitting_out_mean <- merge(model_fitting_out_mean, budburst, by = "tree_id", all.x =T)
model_fitting_out_median <- merge(model_fitting_out_median, budburst, by = "tree_id", all.x =T)
model_fitting_out_sen1 <- merge(model_fitting_out_sen1, budburst, by = "tree_id", all.x =T)

#save results
save(model_fitting_out_all, file = "out/log_function_models/all_values_fitted_models_output.RData")
save(model_fitting_out_mean, file = "out/log_function_models/mean_fitted_models_output.RData")
save(model_fitting_out_median, file = "out/log_function_models/median_fitted_models_output.RData")
save(model_fitting_out_sen1, file = "out/log_function_models/sentinel1_fitted_models_output.RData")

save(models_all, file = "out/log_function_models/all_values_fitted_models.RData")
save(models_mean, file = "out/log_function_models/mean_fitted_models.RData")
save(models_median, file = "out/log_function_models/median_fitted_models.RData")
save(models_sen1, file = "out/log_function_models/sentinel1_fitted_models.RData")



#fit models for NDVI and sigma ratio means over all trees
load("out/all_in_one/sentinel1_daily_sigma_ratio_all_trees.RData")
load("out/all_in_one/aio_daily_ndvi_all_trees_means.RData")
load("out/all_in_one/aio_daily_ndvi_all_trees_medians.RData")

whole_forest_model_fitting <- function(dataset = aio_all_tree_means, return_models = F, sentinel1 = F){
  #perform model fitting and extract SOS, MOS and EOS values from fitted model
  model_fitting_out <- NULL
  models <- list()
  countvar <- 1
  
  aio <- dataset
  
  if(sentinel1 == F){
    #iterate through platforms
    for(platform in unique(aio$platform)){
      
      SOSdoy <- list()
      platform_data_only <- aio[aio$platform == platform,]
      
      ntrees = 50
      SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)
      
      
      #get doys; as not all trees are present in some of the orthomosaics, the doylist is extracted tree-specific
      doylist <- sort(platform_data_only$doy)
      doystart <- head(doylist,1)
      doyend <- tail(doylist,1)
      
      ndvidat <- platform_data_only$ndvi_mean
      
      if(length(ndvidat) == 0){
        print(paste('not data for platform ',platform))
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = NA,
                                                                 no_data = T,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
        
        models[[countvar]] <- list(platform = platform, ndvidat = NA, doylist = NA, model = NA)
        countvar <- countvar+1
      } else {
        
        # see, if model can be calculated; if not, move on
        tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=0.6,b=0.01,c=1,d=0)),error=function(e) e, warning=function(w) w)
        
        if(is(tt,'warning')){
          print(paste('warning at platform ',platform))
          model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                   SOS = NA,
                                                                   MOS = NA,
                                                                   EOS = NA,
                                                                   RSE = NA,
                                                                   no_data = F,
                                                                   warning = T,
                                                                   error = F)) #residualSTDerror
          
          models[[countvar]] <- list(platform = platform, ndvidat = ndvidat, doylist = doylist, model = NA)
          countvar <- countvar+1
        } else if(is(tt,'error')){
          print(paste('error at platform ',platform))
          model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                   SOS = NA,
                                                                   MOS = NA,
                                                                   EOS = NA,
                                                                   RSE = NA,
                                                                   no_data = F,
                                                                   warning = F,
                                                                   error = T)) #residualSTDerror
          
          models[[countvar]] <- list(platform = platform, ndvidat = ndvidat, doylist = doylist, model = NA)
          countvar <- countvar+1
        } else {
          # minpack.lm isneeded for the use of minpack.lm::nlsLM() instead of stats::nls(), as it is more robust to bad starting values (and I couldn't find good ones)
          # it uses the Levenberg-Marquardt (https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) algorithm instead of Gauss-Newton
          fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=0.6,b=0.01,c=1,d=0))
          
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
          predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
          
          #first derivation: slope
          firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2
          
          #second derivation: curvature
          secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2))
          curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
          
          #change of curvature
          ROCcurvature <- c(curvature[2:(doyend-doystart+1)],NA)-curvature 
          
          # diff(ROCcurvature) = amount of change of curvature between two points
          # sign(diff(ROCcurvature)) = positive (1) or negative (-1) change of curvature between two points
          # which(diff(sign(diff(ROCcurvature)))==-2) = point where curvature changes from positive to negative (or the other way around)
          SOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#get DOY of ROC local maximum 1
          MOS <- preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]
          EOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[2]]#get DOY of ROC local maximum 2
          
          model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                   SOS = SOS,
                                                                   MOS = MOS,
                                                                   EOS = EOS,
                                                                   RSE = modsum$sigma,
                                                                   no_data = F,
                                                                   warning = F,
                                                                   error = F)) #residualSTDerror
          
          models[[countvar]] <- list(platform = platform, ndvidat = ndvidat, doylist = doylist, model = fitmodel)
          countvar <- countvar+1
        }
      }
    }
  } else {
    
    SOSdoy <- list()
    ntrees = 50
    SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)
    
    #starting values
    a_start <- 5
    b_start <- 0.5
    c_start <- 120
    d_start <- 5
    
    #get doys; as not all trees are present in some of the orthomosaics, the doylist is extracted tree-specific
    doylist <- sort(aio$doy)
    doystart <- head(doylist,1)
    doyend <- tail(doylist,1)
    
    ndvidat <- aio$sigma_ratio_mean
    
    if(length(ndvidat) == 0){
      model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                               SOS = NA,
                                                               MOS = NA,
                                                               EOS = NA,
                                                               RSE = NA,
                                                               no_data = T,
                                                               warning = F,
                                                               error = F)) #residualSTDerror
      
      models[[countvar]] <- list(platform = "sentinel1", ndvidat = NA, doylist = NA, model = NA)
      countvar <- countvar+1
    } else {
      # see, if model can be calculated; if not, move on
      tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
      
      if(is(tt,'warning')){
        print(paste('warning at platform ',"sentinel1"))
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = NA,
                                                                 no_data = F,
                                                                 warning = T,
                                                                 error = F)) #residualSTDerror
        
        models[[countvar]] <- list(platform = "sentinel1", ndvidat = ndvidat, doylist = doylist, model = NA)
        countvar <- countvar+1
      } else if(is(tt,'error')){
        print(paste('error at platform ',"sentinel1"))
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = NA,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = T)) #residualSTDerror
        
        models[[countvar]] <- list(platform = "sentinel1", ndvidat = ndvidat, doylist = doylist, model = NA)
        countvar <- countvar+1
      } else {
        # minpack.lm isneeded for the use of minpack.lm::nlsLM() instead of stats::nls(), as it is more robust to bad starting values (and I couldn't find good ones)
        # it uses the Levenberg-Marquardt (https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) algorithm instead of Gauss-Newton
        fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start))
        
        #plot(ndvidat ~ doylist)
        #curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)
        
        #get model coefficients
        modsum <- summary(fitmodel)
        a <- modsum$parameters[1]
        b <- modsum$parameters[2]
        c <- modsum$parameters[3]
        d <- modsum$parameters[4]
        
        # predict values of days, that are not present in the images, based on lgistic function "fitmodel"
        preddoylist <- seq(doystart,doyend,1)
        predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
        
        #first derivation: slope
        firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2
        
        #second derivation: curvature
        secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2))
        curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
        
        #change of curvature
        ROCcurvature <- c(curvature[2:(doyend-doystart+1)],NA)-curvature 
        
        # diff(ROCcurvature) = amount of change of curvature between two points
        # sign(diff(ROCcurvature)) = positive (1) or negative (-1) change of curvature between two points
        # which(diff(sign(diff(ROCcurvature)))==-2) = point where curvature changes from positive to negative (or the other way around)
        SOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#get DOY of ROC local maximum 1
        MOS <- preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]
        EOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[2]]#get DOY of ROC local maximum 2
        
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                 SOS = SOS,
                                                                 MOS = MOS,
                                                                 EOS = EOS,
                                                                 RSE = modsum$sigma,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
        
        models[[countvar]] <- list(platform = "sentinel1", ndvidat = ndvidat, doylist = doylist, model = fitmodel)
        countvar <- countvar+1
      }
    }
  }
  model_fitting_out <- model_fitting_out[!duplicated(model_fitting_out[,c(1,2)]),]
  ifelse(return_models == T,return(models),return(model_fitting_out))
}

whole_forest_model_fitting_out_mean <- whole_forest_model_fitting(dataset = aio_all_tree_means, return_models = F, sentinel1 = F)
whole_forest_model_fitting_out_median <- whole_forest_model_fitting(dataset = aio_all_tree_medians, return_models = F, sentinel1 = F)
whole_forest_models_mean <- whole_forest_model_fitting(dataset = aio_all_tree_means, return_models = T, sentinel1 = F)
whole_forest_models_median <- whole_forest_model_fitting(dataset = aio_all_tree_medians, return_models = T, sentinel1 = F)

whole_forest_model_fitting_out_sen1 <- whole_forest_model_fitting(dataset = sen1_all_tree_means, return_models = F, sentinel1 = T)
whole_forest_models_sen1 <- whole_forest_model_fitting(dataset = sen1_all_tree_means, return_models = T, sentinel1 = T)

#save results
save(whole_forest_model_fitting_out_mean, file = "out/log_function_models/whole_forest_mean_fitted_models_output.RData")
save(whole_forest_model_fitting_out_median, file = "out/log_function_models/whole_forest_median_fitted_models_output.RData")
save(whole_forest_model_fitting_out_sen1, file = "out/log_function_models/whole_forest_sentinel1_fitted_models_output.RData")
save(whole_forest_models_mean, file = "out/log_function_models/whole_forest_mean_fitted_models.RData")
save(whole_forest_models_median, file = "out/log_function_models/whole_forest_median_fitted_models.RData")
save(whole_forest_models_sen1, file = "out/log_function_models/whole_forest_sentinel1_fitted_models.RData")
