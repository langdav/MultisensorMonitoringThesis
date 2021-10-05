rm(list = ls())
library(stats);library(rgdal);library(lubridate);library(dplyr);library(RColorBrewer);library(ggplot2);library(minpack.lm) #minpack.lm isneeded for the use of minpack.lm::nlsLM()

#load data
load("data/budburst_data/budburst_phases_def.RData")

#plotting
ggplot(budburst_new[which(budburst_new$tree_id %in% unique(budburst_new$tree_id)[1:5]),]) +
  geom_point(aes(doy, phase_d)) +
  facet_wrap(~tree_id)


model_fitting_out <- NULL
models <- list()
countvar <- 1
SOSdoy <- list()
ntrees = 50
SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)

#starting values
a_start <- 0    #lower asymptote
b_start <- 1   #amplitude of percentage
c_start <- 0.1   #inflection point
d_start <- 0.04   #growth rate

for(phase in colnames(budburst_new)[4:6]){
  for(tree in unique(budburst_new$tree_id)){
    
    #get doys; as not all trees are present in some of the orthomosaics, the doylist is extracted tree-specific
    doylist <- sort(budburst_new$doy[budburst_new$tree_id == tree])
    doystart <- head(doylist,1)
    doyend <- tail(doylist,1)
    
    perc_in_phase <- budburst_new[which(budburst_new$tree_id == tree),phase]
    
    #tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(perc_in_phase ~ b/(1 + exp(-d * (doylist-c))) + a,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
    tt <- tryCatch(fitmodel <- minpack.lm::nlsLM(perc_in_phase ~ b/(1 + exp(c + (doylist*d))) + a,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start)),error=function(e) e, warning=function(w) w)
    
    if(is(tt,'warning')){
      print(paste('warning at tree ',tree,'in',phase))
      model_fitting_out <- rbind(model_fitting_out, data.frame(phase = phase,
                                                               tree_id = tree,
                                                               SOS = NA,
                                                               MOS = NA,
                                                               EOS = NA,
                                                               RSE = NA,
                                                               no_data = F,
                                                               warning = T,
                                                               error = F)) #residualSTDerror
      
      models[[countvar]] <- list(phase = phase, 
                                 tree_id = tree, 
                                 perc_in_phase = perc_in_phase, 
                                 doylist = doylist, 
                                 model = NA)
      countvar <- countvar+1
    } else if(is(tt,'error')){
      print(paste('error at tree ',tree,'in',phase))
      model_fitting_out <- rbind(model_fitting_out, data.frame(phase = phase,
                                                               tree_id = tree,
                                                               SOS = NA,
                                                               MOS = NA,
                                                               EOS = NA,
                                                               RSE = NA,
                                                               no_data = F,
                                                               warning = F,
                                                               error = T)) #residualSTDerror
      
      models[[countvar]] <- list(phase = phase,
                                 tree_id = tree,
                                 perc_in_phase = perc_in_phase,
                                 ndvi_sd = NA,
                                 doylist = doylist,
                                 model = NA)
      countvar <- countvar+1
    } else {
      # minpack.lm isneeded for the use of minpack.lm::nlsLM() instead of stats::nls(), as it is more robust to bad starting values (and I couldn't find good ones)
      # it uses the Levenberg-Marquardt (https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm) algorithm instead of Gauss-Newton
      #fitmodel <- minpack.lm::nlsLM(perc_in_phase ~ b/(1 + exp(-d * (doylist-c))) + a,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start))
      fitmodel <- minpack.lm::nlsLM(perc_in_phase ~ b/(1 + exp(c + (doylist*d))) + a,control=nls.control(warnOnly=TRUE,maxiter = 100),start=list(a=a_start,b=b_start,c=c_start,d=d_start))
      
      # predict values of days, that are not present in the images, based on lgistic function "fitmodel"
      preddoylist <- seq(doystart,doyend,1)
      predNDVImod <- predict(fitmodel,data.frame(doylist=preddoylist))
      
      # ggplot() +
      #   geom_point(aes(doylist, perc_in_phase)) +
      #   geom_line(aes(preddoylist,predNDVImod))
      
      #get model coefficients
      modsum <- summary(fitmodel)
      a <- modsum$parameters[1]
      b <- modsum$parameters[2]
      c <- modsum$parameters[3]
      d <- modsum$parameters[4]
      
      #first derivation: slope
      #firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2 # old model
      firstDerivNDVImod <- -1*((b*d*exp(c + (preddoylist*d)))/(exp(c + (preddoylist*d))+1)^2) # new model
      
      #second derivation: curvature
      # secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2)) #old model
      secondDerivNDVImod <- ((2*b*d^2*exp(2*c + (2*preddoylist*d)))/((exp(c + (preddoylist*d))+1)^3))-((b*d^2*exp(c + (preddoylist*d)))/((exp(c + (preddoylist*d))+1)^2)) #new model
      curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
      
      #curvature change rate (= first derivation)
      # firstDerivCurv <- (b*(d^3)*exp(d*(preddoylist-c))*(exp(6*d*(preddoylist-c))+(-2*(b^2)*(d^2)-9)*exp(4*d*(preddoylist-c))+(2*(b^2)*(d^2)-16)*exp(3*d*(preddoylist-c))+(-2*(b^2)*(d^2)-9)*exp(2*d*(preddoylist-c))+1)) / (((((b^2)*(d^2)*exp(-2*d*(preddoylist-c)))/((exp(-1*d*(preddoylist-c))+1)^4)+1)^2.5)*((exp(d*(preddoylist-c))+1)^8)) #old model
      firstDerivCurv <- -1*(b*(d^3)*exp(c+(preddoylist*d))*(exp(6*c+(6*preddoylist*d))+(-2*(b^2)*exp(4*c)*(d^2)-9*exp(4*c))*exp(4*d*preddoylist)+(2*(b^2)*exp(3*c)*(d^2)-16*exp(3*c))*exp(3*d*preddoylist)+(-2*(b^2)*exp(2*c)*(d^2)-9*exp(2*c))*exp(2*d*preddoylist)+1) / (((exp(c+(preddoylist*d))+1)^8)*((((((b^2)*(d^2)*exp(2*(c+(preddoylist*d)))) / ((exp(c+(preddoylist*d))+1)^4))) + 1)^2.5))) #new model
      
      #change of curvature
      ROCcurvature <- c(curvature[2:length(curvature)],NA)-curvature
      
      # old way of finding SOS
      # SOS <- sort(preddoylist[kit::topn(ROCcurvature,2)])[1]
      # MOS <- preddoylist[which.min(ROCcurvature)]
      # EOS <- sort(preddoylist[kit::topn(ROCcurvature,2)])[2]
      # 
      
      if(all(is.na(firstDerivCurv)==T)){
        model_fitting_out <- rbind(model_fitting_out, data.frame(phase = phase,
                                                                 tree_id = tree,
                                                                 SOS = NA,
                                                                 MOS = NA,
                                                                 EOS = NA,
                                                                 RSE = modsum$sigma,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
      } else {
        # SOS = first max. of firstDerivCurv; MOS = min. of firstDerivCurv; EOS = second max. of firstDerivCurv; 
        SOS <- sort(preddoylist[kit::topn(firstDerivCurv,2)])[1]
        MOS <- preddoylist[which.min(firstDerivCurv)]
        EOS <- sort(preddoylist[kit::topn(firstDerivCurv,2)])[2]
        
        
        model_fitting_out <- rbind(model_fitting_out, data.frame(phase = phase,
                                                                 tree_id = tree,
                                                                 SOS = SOS,
                                                                 MOS = MOS,
                                                                 EOS = EOS,
                                                                 RSE = modsum$sigma,
                                                                 no_data = F,
                                                                 warning = F,
                                                                 error = F)) #residualSTDerror
      }
      
      models[[countvar]] <- list(phase = phase,
                                 tree_id = tree,
                                 perc_in_phase = perc_in_phase,
                                 ndvi_sd = NA,
                                 doylist = doylist,
                                 model = fitmodel)
      
      countvar <- countvar+1
    }
  }
}

budburst_model_fitting_out <- model_fitting_out;rm(model_fitting_out)
budbust_models <- models;rm(models)

#save models
save(budburst_model_fitting_out, file = "out/log_function_models/budburst_fitted_models_output.RData")
save(budbust_models, file = "out/log_function_models/budburst_fitted_models.RData")

#plotting
plot_SOS <- function(phase = "phase_d", tree = 10){
  #x- and y-axis limits
  xlim_min <- min(unique(budburst_new$doy))
  xlim_max <- max(unique(budburst_new$doy))
  ylim_min <- -0.1
  ylim_max <- 1.1
  
  #select platform; one of sentinel2, orthomosaic, treetalker, planetscope
  phase_sel <- phase
  
  #select tree; between 1 and 50
  tree_sel <- unique(budburst_new$tree_id)[tree]
  
  entry <- NULL
  for(i in 1:length(budbust_models)){
    if(budbust_models[[i]][["tree_id"]] == tree_sel & budbust_models[[i]][["phase"]] == phase_sel){
      entry <- c(entry, i)
    }
  }
  model_sel <- budbust_models[[entry]]
  data_sel <- budburst_model_fitting_out[which(budburst_model_fitting_out$phase == phase_sel & budburst_model_fitting_out$tree_id == tree_sel),]
  
  perc_in_phase <- model_sel$perc_in_phase
  fitmodel <- model_sel$model
  if(!any(is.na(fitmodel))){
    doylist <- model_sel$doylist
    phase <- model_sel$phase
    tree <- model_sel$tree
    preddoylist <- seq(head(doylist, 1),tail(doylist,1),1)
    predNDVImod <- predict(fitmodel,data.frame(doylist=preddoylist))
    
    if(!is.na(data_sel$SOS)){
      plot_out <- ggplot() +
        geom_point(aes(doylist, perc_in_phase)) +
        geom_line(aes(preddoylist,predNDVImod)) +
        geom_vline(xintercept = data_sel$SOS, linetype = "dotdash", color = "blue", size = .8) +
        geom_text(aes(x = data_sel$SOS,
                      y = min(perc_in_phase),
                      label = "predicted",
                      vjust = 0,
                      angle = 90)) +
        xlim(xlim_min, xlim_max) +
        ylim(ylim_min, ylim_max) +
        theme_light() +
        xlab("Day of Year") +
        ylab(paste0("% of tree in",phase)) +
        ggtitle(paste0(phase_sel,";",tree))
    } else {
      plot_out <- ggplot() +
        geom_point(aes(doylist, perc_in_phase)) +
        geom_line(aes(preddoylist,predNDVImod)) +
        xlim(xlim_min, xlim_max) +
        ylim(ylim_min, ylim_max) +
        theme_light() +
        xlab("Day of Year") +
        ylab(paste0("% of tree in",phase)) +
        ggtitle(paste0(phase_sel,";",tree))
    }
  } else {
    doylist <- model_sel$doylist
    phase <- model_sel$phase
    tree <- model_sel$tree
    
    plot_out <- ggplot() +
      geom_point(aes(doylist, perc_in_phase)) +
      xlim(xlim_min, xlim_max) +
      ylim(ylim_min, ylim_max) +
      theme_light() +
      xlab("Day of Year") +
      ylab(paste0("% of tree in",phase)) +
      ggtitle(paste0(phase_sel,";",tree))
  }
  return(plot_out)
}

#plot single tree in single phase
plot_SOS(phase = "phase_d", tree = 10)

#plot all trees and all phases
for(phase in colnames(budburst_new)[4:6]){
  all_trees <- lapply(1:25, function(x) plot_SOS(phase = phase, tree = x))
  all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
  #all_trees_panel
  ggsave(paste0("out/model_plots/budburst_obs_",phase,"_1_25",".png"),all_trees_panel,width = 20,height=12)
  
  all_trees <- lapply(26:50, function(x) plot_SOS(phase = phase, tree = x))
  all_trees_panel <- gridExtra::marrangeGrob(all_trees, nrow = 5, ncol = 5)
  #all_trees_panel
  ggsave(paste0("out/model_plots/budburst_obs_",phase,"_25_50",".png"),all_trees_panel,width = 20,height=12)
}
