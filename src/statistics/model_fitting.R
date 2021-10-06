#author: David Langenohl; based on the work of Dominic Fawcett
#last modified: 19.08.2021
#description: script to extract NDVI time series per crown and fit logistic function to retrieve phenological start and end of budburst phase
#NOTE: no models fitted to Sentinel-2 data, as too little data available (3 days)
#NOTE: model variables: a = amplitude of increase; b = growth rate; c = value of t at inflection point; d = lower asymptote

rm(list = ls())
library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);library(kit);library(minpack.lm) #minpack.lm isneeded for the use of minpack.lm::nlsLM()

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
      # if(platform == "planetscope"){
      #   a_start <- 0.7    #inflection point
      #   b_start <- 0.01   #growth rate
      #   c_start <- 1      #amplitude of increase
      #   d_start <- 0.5    #lower asymptote
      # } else if(platform == "sentinel2"){
      #   a_start <- 0.7
      #   b_start <- 0.01
      #   c_start <- 1
      #   d_start <- 0.5
      # } else if(platform == "treetalker"){
      #   a_start <- 0.7
      #   b_start <- 0.01
      #   c_start <- 1
      #   d_start <- 0
      # } else if(platform == "orthomosaic"){
      #   a_start <- 120
      #   b_start <- 0.001
      #   c_start <- 1
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
          sd <- platform_data_only$ndvi_sd[which(platform_data_only$tree_id == tree)]
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
          
          models[[countvar]] <- list(platform = platform,
                                     tree_id = tree,
                                     ndvidat = NA,
                                     ndvi_sd <- if(mean == T & median == F){sd} else {NA},
                                     doylist = NA,
                                     model = NA)
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
            
            models[[countvar]] <- list(platform = platform, 
                                       tree_id = tree,
                                       ndvidat = ndvidat,
                                       ndvi_sd <- if(mean == T & median == F){sd} else {NA},
                                       doylist = doylist, 
                                       model = NA)
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
            
            models[[countvar]] <- list(platform = platform, 
                                       tree_id = tree, 
                                       ndvidat = ndvidat,
                                       ndvi_sd <- if(mean == T & median == F){sd} else {NA},
                                       doylist = doylist, 
                                       model = NA)
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
            predNDVImod <- predict(fitmodel,data.frame(doylist=preddoylist))
            
            #first derivation: slope
            # firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2 # old model
            firstDerivNDVImod <- -1*((b*d*exp(c + (preddoylist*d)))/(exp(c + (preddoylist*d))+1)^2) # new model
            
            #second derivation: curvature
            # secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2)) #old model
            secondDerivNDVImod <- ((2*b*d^2*exp(2*c + (2*preddoylist*d)))/((exp(c + (preddoylist*d))+1)^3))-((b*d^2*exp(c + (preddoylist*d)))/((exp(c + (preddoylist*d))+1)^2)) #new model
            curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
            
            #curvature change rate (= first derivation)
            #firstDerivCurv <- (b*(d^3)*exp(d*(preddoylist-c))*(exp(6*d*(preddoylist-c))+(-2*(b^2)*(d^2)-9)*exp(4*d*(preddoylist-c))+(2*(b^2)*(d^2)-16)*exp(3*d*(preddoylist-c))+(-2*(b^2)*(d^2)-9)*exp(2*d*(preddoylist-c))+1)) / (((((b^2)*(d^2)*exp(-2*d*(preddoylist-c)))/((exp(-1*d*(preddoylist-c))+1)^4)+1)^2.5)*((exp(d*(preddoylist-c))+1)^8)) #old model
            firstDerivCurv <- -1*(b*(d^3)*exp(c+(preddoylist*d))*(exp(6*c+(6*preddoylist*d))+(-2*(b^2)*exp(4*c)*(d^2)-9*exp(4*c))*exp(4*d*preddoylist)+(2*(b^2)*exp(3*c)*(d^2)-16*exp(3*c))*exp(3*d*preddoylist)+(-2*(b^2)*exp(2*c)*(d^2)-9*exp(2*c))*exp(2*d*preddoylist)+1) / (((exp(c+(preddoylist*d))+1)^8)*((((((b^2)*(d^2)*exp(2*(c+(preddoylist*d)))) / ((exp(c+(preddoylist*d))+1)^4))) + 1)^2.5))) #new model
            
            # ggplot() + 
            #   geom_line(aes(preddoylist, firstDerivCurv)) + 
            #   scale_x_continuous(breaks = round(seq(min(doylist), max(doylist), by = 1),1))
            
            #change of curvature
            ROCcurvature <- c(curvature[2:length(curvature)],NA)-curvature
            
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
            
            
            
            # old way of finding SOS
            # SOS <- sort(preddoylist[kit::topn(ROCcurvature,2)])[1]
            # MOS <- preddoylist[which.min(ROCcurvature)]
            # EOS <- sort(preddoylist[kit::topn(ROCcurvature,2)])[2]
            
            if(all(is.na(firstDerivCurv)==T)){
              model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
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
              
              model_fitting_out <- rbind(model_fitting_out, data.frame(platform = platform,
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
              #   geom_line(aes(preddoylist,firstDerivCurv*1000+0.7), linetype = "dashed") +
              #   geom_vline(xintercept = SOS, linetype = "longdash", colour = "red") +
              #   geom_text(aes(x = SOS-1,
              #                 y = 0.4,
              #                 label = "budburst",
              #                 vjust = 0,
              #                 angle = 90)) +
              #   theme(panel.grid = element_blank(),
              #         axis.title.y = element_blank(),
              #         axis.text.y = element_blank(),
              #         axis.ticks.y = element_blank(),
              #         panel.background = element_blank()) +
              #   xlab("DOY")
              # ggsave(paste0("tmp/example_curvature_derivative.png"),width = 5,height=4)
            }
            
            models[[countvar]] <- list(platform = platform, 
                                       tree_id = tree,
                                       ndvidat = ndvidat,
                                       ndvi_sd <- if(mean == T & median == F){sd} else {NA},
                                       doylist = doylist,
                                       model = fitmodel)
            
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
    # a_start <- 6    #inflection point
    # b_start <- 0.01   #growth rate
    # c_start <- 10      #amplitude of increase
    # d_start <- 5    #lower asymptote
    
    a_start <- 5    #lower asymptote
    b_start <- 10   #amplitude of backscatter
    c_start <- 6    #inflection point
    d_start <- 0.01 #growth rate
    
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
        
        models[[countvar]] <- list(platform = "sentinel1", 
                                   tree_id = tree, 
                                   ndvidat = NA, 
                                   ndvi_sd = NA,
                                   doylist = NA, 
                                   model = NA)
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
          
          models[[countvar]] <- list(platform = "sentinel1", 
                                     tree_id = tree,
                                     ndvidat = ndvidat,
                                     ndvi_sd = NA,
                                     doylist = doylist,
                                     model = NA)
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
          predNDVImod <- predict(fitmodel,data.frame(doylist=preddoylist))
          
          #first derivation: slope
          #firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2 # old model
          firstDerivNDVImod <- -1*((b*d*exp(c + (preddoylist*d)))/(exp(c + (preddoylist*d))+1)^2) # new model
          
          #second derivation: curvature
          #secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2)) #old model
          secondDerivNDVImod <- ((2*b*d^2*exp(2*c + (2*preddoylist*d)))/((exp(c + (preddoylist*d))+1)^3))-((b*d^2*exp(c + (preddoylist*d)))/((exp(c + (preddoylist*d))+1)^2)) #new model
          curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
          
          #curvature change rate (= first derivation)
          #firstDerivCurv <- (b*(d^3)*exp(d*(preddoylist-c))*(exp(6*d*(preddoylist-c))+(-2*(b^2)*(d^2)-9)*exp(4*d*(preddoylist-c))+(2*(b^2)*(d^2)-16)*exp(3*d*(preddoylist-c))+(-2*(b^2)*(d^2)-9)*exp(2*d*(preddoylist-c))+1)) / (((((b^2)*(d^2)*exp(-2*d*(preddoylist-c)))/((exp(-1*d*(preddoylist-c))+1)^4)+1)^2.5)*((exp(d*(preddoylist-c))+1)^8)) #old model
          firstDerivCurv <- -1*(b*(d^3)*exp(c+(preddoylist*d))*(exp(6*c+(6*preddoylist*d))+(-2*(b^2)*exp(4*c)*(d^2)-9*exp(4*c))*exp(4*d*preddoylist)+(2*(b^2)*exp(3*c)*(d^2)-16*exp(3*c))*exp(3*d*preddoylist)+(-2*(b^2)*exp(2*c)*(d^2)-9*exp(2*c))*exp(2*d*preddoylist)+1) / (((exp(c+(preddoylist*d))+1)^8)*((((((b^2)*(d^2)*exp(2*(c+(preddoylist*d)))) / ((exp(c+(preddoylist*d))+1)^4))) + 1)^2.5))) #new model
          
          #change of curvature
          ROCcurvature <- c(curvature[2:length(curvature)],NA)-curvature
          
          # old way of finding SOS
          # SOS <- sort(preddoylist[kit::topn(ROCcurvature,2)])[1]
          # MOS <- preddoylist[which.min(ROCcurvature)]
          # EOS <- sort(preddoylist[kit::topn(ROCcurvature,2)])[2]
          
          if(all(is.na(firstDerivCurv)==T)){
            model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
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
            
            model_fitting_out <- rbind(model_fitting_out, data.frame(platform = "sentinel1",
                                                                     tree_id = tree,
                                                                     SOS = SOS,
                                                                     MOS = MOS,
                                                                     EOS = EOS,
                                                                     RSE = modsum$sigma,
                                                                     no_data = F,
                                                                     warning = F,
                                                                     error = F)) #residualSTDerror
          }
          
          models[[countvar]] <- list(platform = "sentinel1", 
                                     tree_id = tree,
                                     ndvidat = ndvidat,
                                     ndvi_sd = NA,
                                     doylist = doylist,
                                     model = fitmodel)
          
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
i <- 165
ndvidat <- models_mean[[i]]$ndvidat;doylist <- models_mean[[i]]$doylist;fitmodel <- models_mean[[i]]$model;plot(ndvidat ~ doylist, main = models_mean[[i]]$platform);curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)  

#sentinel-1
i <- 4
ndvidat <- models_sen1[[i]]$ndvidat;doylist <- models_sen1[[i]]$doylist;fitmodel <- models_sen1[[i]]$model;plot(ndvidat ~ doylist, main = models_sen1[[i]]$platform);curve(predict(fitmodel, newdata = data.frame(doylist = x)), add = TRUE)  



#for each tree, add doy of manually observed budburst
# budburst <- read.csv("data/budburst_data/budburst_long.csv")
# budburst$budburst_obervation_doy <- yday(budburst$date)
# budburst <- budburst[which(budburst$budburst==T),c("tree_id","budburst_obervation_doy","budburst_perc")]
# budburst <- budburst[!duplicated(budburst$tree_id),]
# 
# model_fitting_out_all <- merge(model_fitting_out_all, budburst, by = "tree_id", all.x =T)
# model_fitting_out_mean <- merge(model_fitting_out_mean, budburst, by = "tree_id", all.x =T)
# model_fitting_out_median <- merge(model_fitting_out_median, budburst, by = "tree_id", all.x =T)
# model_fitting_out_sen1 <- merge(model_fitting_out_sen1, budburst, by = "tree_id", all.x =T)

#add predicted budburst from observations
load("out/log_function_models/budburst_fitted_models_output.RData")
obs_budburst <- merge(budburst_model_fitting_out[which(budburst_model_fitting_out$phase == "phase_d"),c(2,3)],
                      budburst_model_fitting_out[which(budburst_model_fitting_out$phase == "phase_e"),c(2,3)],by="tree_id")
obs_budburst <- merge(obs_budburst, budburst_model_fitting_out[which(budburst_model_fitting_out$phase == "phase_f"),c(2,3)],by="tree_id")
colnames(obs_budburst) <- c("tree_id", "SOS_phase_d", "SOS_phase_e", "SOS_phase_f")

model_fitting_out_all <- merge(model_fitting_out_all, obs_budburst, by = "tree_id", all.x =T)
model_fitting_out_mean <- merge(model_fitting_out_mean, obs_budburst, by = "tree_id", all.x =T)
model_fitting_out_median <- merge(model_fitting_out_median, obs_budburst, by = "tree_id", all.x =T)
model_fitting_out_sen1 <- merge(model_fitting_out_sen1, obs_budburst, by = "tree_id", all.x =T)

#save results
save(model_fitting_out_all, file = "out/log_function_models/all_values_fitted_models_output.RData")
save(model_fitting_out_mean, file = "out/log_function_models/mean_fitted_models_output.RData")
save(model_fitting_out_median, file = "out/log_function_models/median_fitted_models_output.RData")
save(model_fitting_out_sen1, file = "out/log_function_models/sentinel1_fitted_models_output.RData")

save(models_all, file = "out/log_function_models/all_values_fitted_models.RData")
save(models_mean, file = "out/log_function_models/mean_fitted_models.RData")
save(models_median, file = "out/log_function_models/median_fitted_models.RData")
save(models_sen1, file = "out/log_function_models/sentinel1_fitted_models.RData")
