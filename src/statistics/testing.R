library(stats)
stats::nls()


library(rgdal)
library(lubridate)
library(RColorBrewer)
library(ggplot2)
library(raster)

##########################
# inital parameters for beech for successful fitting

# Orthomosaic
load("out/orthomosaic/ndvi_long_format_phenoclasses_orthomosaic.RData") #only load if necessary, as resulting dataframe is huge and contains 8553222 rows
orthomosaic <- ndvi_long_pheno_orthomosaic; rm(ndvi_long_pheno_orthomosaic)
ortho <- orthomosaic[which(orthomosaic$tree_id == "mof_cst_00001"),]

load("out/all_in_one/aio_daily_ndvi_per_tree_means.RData")
aio <- aio_daily_ndvi_means;rm(aio_daily_ndvi_means)
ortho <- aio[aio$platform == "planetscope",]
#perform model fitting and extract SOS, MOS and EOS values from fitted model

#get doys
doylist <- yday(unique(ortho$date[ortho$tree_id == "mof_cst_00001"]))
doystart <- sort(doylist)[1]
doyend <- sort(doylist)[length(doylist)]


#iterate through classes
for(i in 1:18){
  #i=1#debug
  
  SOSdoy <- list()
  classNDVIvals <- ortho
  
  ntrees <- length(unique(ortho$tree_id))
  doylistcurrent <- doylist
  
  SOSdoymat <- matrix(NA,nrow=1,ncol=ntrees)
  
  #iterate through individuals per class
  for(j in c(1:ntrees)){
    j <- 1
    tree <- unique(ortho$tree_id)[j]
    
    #ndvidat <- unlist(classNDVIvals[which(classNDVIvals$tree_id == tree),])
    ndvidat <- classNDVIvals$ndvi_mean[which(classNDVIvals$tree_id == tree)]
    
    
    #inital parameters for beech for successful fitting
    fitmodel <- nls(ndvidat ~ a/(1 + exp(-b * (doylistcurrent-c)) ) + d,control=nls.control(warnOnly=TRUE,maxiter = 50),start=list(a=0.6,b=0.01,c=1,d=0.3))
    
    #get model coefficients
    modsum <- summary(fitmodel)
    a <- modsum$parameters[1]
    b <- modsum$parameters[2]
    c <- modsum$parameters[3]
    d <- modsum$parameters[4]
    preddoylist <- seq(doystart,doyend,1)
    predNDVImod <- predict(fitmodel,data.frame(doylistcurrent=preddoylist))
    
    
    #slope
    firstDerivNDVImod <- (a*b*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2
    #curvature
    secondDerivNDVImod <- a*(((2*b^2*exp(-2*b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^3)-((b^2*exp(-b * (preddoylist-c)))/(1+exp((-b * (preddoylist-c))))^2))
    
    curvature <- secondDerivNDVImod/((1+(firstDerivNDVImod)^2)^1.5)
    
    ROCcurvature <- c(curvature[2:(doyend-doystart+1)],NA)-curvature
    
    if(length((which(diff(sign(diff(ROCcurvature)))==-2)+1))<2){ #if SOS can't be retrieved, ignore MOS and EOS
      print(paste('SOS fail at class', i,', individual ',j))
      SOS <- NA
      MOS <- NA
      EOS <- NA
    }else{ 
      SOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[1]]#get DOY of ROC local maximum 1
      MOS <- preddoylist[firstDerivNDVImod==max(firstDerivNDVImod,na.rm=TRUE)]
      EOS <- preddoylist[(which(diff(sign(diff(ROCcurvature)))==-2)+1)[2]]#get DOY of ROC local maximum 2
    }
    
    SOSdoylist[indcount,1] <- SOS
    
    SOSdoylist[indcount,2] <- MOS
    
    SOSdoylist[indcount,3] <- EOS
    SOSdoylist[indcount,4] <- i
    
    residualSTDerror[indcount,1] <- modsum$sigma
    residualSTDerror[indcount,2] <- i
    indcount <- indcount+1
    
    
  }
  
}

#make data frame with extracted dates per individual and class
SOSdf <- data.frame(SOS=SOSdoylist[,1],MOS=SOSdoylist[,2],EOS=SOSdoylist[,3],class=SOSdoylist[,4],RSE=residualSTDerror[,1])




