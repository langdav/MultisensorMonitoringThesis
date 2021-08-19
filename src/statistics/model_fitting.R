#author: David Langenohl; based on the work of Dominic Fawcett
#last modified: 19.08.2021
#description: script to extract NDVI time series per crown and fit logistic function to retrieve phenological start and end of budburst phase
#NOTE: no models fitted to Sentinel-2 data, as too little data available (3 days)

rm(list = ls())
library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);library(minpack.lm) #minpack.lm isneeded for the use of minpack.lm::nlsLM()

# load data
load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
aio <- aio[-which(aio$date > as.Date("2021-07-01")),] #limit to data before "2021-07-01"
#aio <- aio[-which(aio$platform == "sentinel2"),] #remove Sentinel-2 data, as there are only 3 dates available which are not enough to fit a proper model

#perform model fitting and extract SOS, MOS and EOS values from fitted model
model_fitting_out <- NULL
models <- list()
countvar <- 1

#iterate through platforms
for(platform in unique(aio$platform)){
  
  SOSdoy <- list()
  platform_data_only <- aio[aio$platform == platform,]
  
  ntrees = 50
  SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)
  
  #iterate through all trees per platform
  for(tree in unique(aio$tree_id)){
    
    #get doys; as not all trees are present in some of the orthomosaics, the doylist is extracted tree-specific
    doylist <- yday(unique(platform_data_only$date[platform_data_only$tree_id == tree]))
    doystart <- head(doylist,1)
    doyend <- tail(doylist,1)
    
    #ndvidat <- unlist(classNDVIvals[which(classNDVIvals$tree_id == tree),])
    ndvidat <- platform_data_only$ndvi_mean[which(platform_data_only$tree_id == tree)]
    
    if(length(ndvidat) == 0){
      print(paste('not data for platform ',platform, ', tree ',tree))
      model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                               tree = tree,
                                                               SOS = NA,
                                                               MOS = NA,
                                                               EOS = NA,
                                                               RSE = NA,
                                                               no_data = T,
                                                               warning = F,
                                                               error = F)) #residualSTDerror
      
      models[[countvar]] <- list(platform = platform, tree = tree, ndvidat = NA, doylist = NA, model = NA)
      countvar <- countvar+1
    } else {
      
      # see, if model can be calculated; if not, move on
      tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(ndvidat ~ a/(1 + exp(-b * (doylist-c))) + d,control=nls.control(warnOnly=TRUE,maxiter = 500),start=list(a=0.6,b=0.01,c=1,d=0)),error=function(e) e, warning=function(w) w)
      
      if(is(tt,'warning')){
        print(paste('warning at platform ',platform, ', tree ',tree))
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                 tree = tree,
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = NA,
                                                                 no_data = F,
                                                                 warning = T,
                                                                 error = F)) #residualSTDerror
        
        models[[countvar]] <- list(platform = platform, tree = tree, ndvidat = ndvidat, doylist = doylist, model = NA)
        countvar <- countvar+1
      } else if(is(tt,'error')){
        print(paste('error at platform ',platform, ', tree ',tree))
        model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
                                                                 tree = tree,
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = NA,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = T)) #residualSTDerror
        
        models[[countvar]] <- list(platform = platform, tree = tree, ndvidat = ndvidat, doylist = doylist, model = NA)
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
                                                                 tree = tree,
                                                                 SOS = SOS,
                                                                 MOS = MOS,
                                                                 EOS = EOS,
                                                                 RSE = modsum$sigma,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
        
        models[[countvar]] <- list(platform = platform, tree = tree, ndvidat = ndvidat, doylist = doylist, model = fitmodel)
        countvar <- countvar+1
      }
    }
  }
}

save(model_fitting_out, file = "out/log_function_models/fitted_models_output.RData")
save(models, file = "out/log_function_models/fitted_models.RData")
