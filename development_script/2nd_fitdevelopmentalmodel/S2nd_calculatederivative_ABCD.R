## This script is to generate derivative and posterior derivative values from GAMM models.
## Age will be sampled from 8.9 to 13.5, with a total of 1000 data points. 
## Covariates will be set as median or mode.
## The predicted derivatives will be generated.
## The data will be used to draw developmental trajectories and analyse the developmental alignment.
library(tidyverse)
library(R.matlab)
library(psych)
library(gratia)
library(mgcv)
library(parallel)

rm(list = ls())
# set resolution
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
# set path
wdpath <- getwd()
if (str_detect(wdpath, "Users")){
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
}else{
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_ABCD'
}

#### load data
CVthr = 75
gamresultsum<-readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
SCrank <- matrix(NA, ds.resolution, ds.resolution)
for (x in 1:ds.resolution){
  for (y in 1:ds.resolution){
    SCrank[x,y] = x^2 + y^2
  }
}
SCrank <- SCrank[lower.tri(SCrank, diag = T)]
SCrank.df <- data.frame(parcel=paste0("SC.", 1:elementnum, "_h"), SCrank=SCrank)
#### source function
source(paste0(functionFolder, '/gamderivatives.R'))

## derivative
#######################################
derivative.sum<-mclapply(1:nrow(gamresultsum), function(x){
  SClabel.tmp <- gamresultsum$parcel[x]
  modobj<-gammodelsum[[x]]
  draws<-1
  increments<-1000
  derivdata<- gam.derivatives(modobj, "age",draws, increments, return_posterior_derivatives = FALSE)
  derivdata$label_ID<-gamresultsum$parcel[x]
  meanSC <- mean(modobj$gam$model[,SClabel.tmp],na.rm=T)
  derivdata$meanSC<-meanSC
  return(derivdata)
}, mc.cores = 30)

derivative.df<-do.call(rbind, lapply(derivative.sum, function(x) data.frame(x)))
saveRDS(derivative.df, paste0(resultFolder, '/derivative.df', elementnum, '_CV', CVthr,'.rds'))

## posterior derivative
############################
derivative.posterior.sum<-mclapply(1:nrow(gamresultsum), function(x){
  modobj<-gammodelsum[[x]]
  draws<-1000
  increments<-1000
  SClabel<-gamresultsum$parcel[x]
  derivdata<- gam.derivatives(modobj, "age",draws, increments, return_posterior_derivatives = TRUE)
  derivdata$SCrank<-SCrank.df$SCrank[SCrank.df$parcel==SClabel]
  meanSC <- mean(modobj$gam$model[,SClabel],na.rm=T)
  derivdata$meanSC<-meanSC
  maxSC <- max(modobj$gam$model[,SClabel],na.rm=T)
  derivdata$maxSC<-maxSC
  return(derivdata)
}, mc.cores = 30)
saveRDS(derivative.posterior.sum, paste0(resultFolder, '/derivative.posterior.df.SA', ds.resolution,'_CV', CVthr,'.rds'))
