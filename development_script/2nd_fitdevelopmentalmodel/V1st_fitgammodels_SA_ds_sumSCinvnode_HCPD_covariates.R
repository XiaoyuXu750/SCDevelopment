## This script is to fit gam models for each edge.
## Sensitivity analysis for additional covariates.
rm(list=ls())
library(mgcv)
library(parallel)
library(tidyverse)
wdpath <- getwd()
# set resolution
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
addcovariate <- "income.adj" # income.adj / ICV
# set path
if (str_detect(wdpath, "Users")){
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
  functionFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
}else{
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_HCPD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
}
# set which consistency threshold should be used
CVthr=75
# load data
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata_SA", ds.resolution,"_CV", CVthr,"_sumSCinvnode.sum.msmtcsd.combatage.rds"))
nrow(SCdata.sum.merge)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
SCdata.sum.merge$totalstrength <- rowMeans(SCdata.sum.merge[,str_detect(names(SCdata.sum.merge), "SC.") & str_detect(names(SCdata.sum.merge), "_h")])
summary(SCdata.sum.merge$totalstrength)
summary(SCdata.sum.merge$ICV)
# source function
source(paste0(functionFolder, '/gamsmooth.R'))
source(paste0(functionFolder, '/plotdata_generate.R'))
detectCores()
## calculate gam results
covariates<-paste0("sex+mean_fd+", addcovariate)
dataname<-"SCdata.sum.merge"
smooth_var<-"age"
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-grep("SC.", names(SCdata.sum.merge), value=T)[x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
    gamresult<-as.data.frame(gamresult)
    return(gamresult)
  }, mc.cores = 50)
  
  gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
  gamresultsum.df[,c(2:18)]<-lapply(gamresultsum.df[,c(2:18)], as.numeric)
  class(gamresultsum.df$gam.smooth.pvalue)
  class(gamresultsum.df$partialRsq)
  # add fdr p values
  gamresultsum.df$pfdr<- p.adjust(gamresultsum.df$anova.smooth.pvalue, method = "fdr")
  gamresultsum.df$sig<-(gamresultsum.df$pfdr<0.05)
  summary(gamresultsum.df)
  saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_validat_', addcovariate, '.rds'))
}else{gamresultsum.df <- readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_validat_', addcovariate, '.rds'))}
## calculate gam models
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-names(SCdata.sum.merge)[1+x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    return(gamresult)
  }, mc.cores = 50)
  saveRDS(resultsum, paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_validat_', addcovariate, '.rds'))
}else{resultsum <- readRDS(paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_validat_', addcovariate, '.rds'))}

## plot data raw
#### generate function data
if (str_detect(wdpath, "cuizaixu_lab")){
  gammodelsum <- resultsum
  plotdatasum<-mclapply(1:elementnum, function(x){
    modobj<-gammodelsum[[x]]
    plotdata<- plotdata_generate(modobj, "age")
    plotdata$SC_label <- names(plotdata)[15]
    plotdata[,15] <- NULL
    return(plotdata)
  }, mc.cores = 50)
  plotdatasum.df <- do.call(rbind, lapply(plotdatasum, function(x) data.frame(x)))
  saveRDS(plotdatasum.df, paste0(interfileFolder, '/plotdatasum.df_SA', ds.resolution,'_sumSCinvnode_CV', CVthr,'_validat_', addcovariate, '.rds'))
}else{plotdatasum.df <- readRDS(paste0(interfileFolder, '/plotdatasum.df_SA', ds.resolution,'_sumSCinvnode_CV', CVthr,'_validat_', addcovariate, '.rds'))}

# To avoid the influence of averaged weight on derivative analyses, we divided SC strength of each edges by their
# weight at age of 8.
SCdata.diw <- SCdata.sum.merge
for (i in 1:elementnum){
  SClabel <- names(SCdata.sum.merge)[i+1]
  plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label==SClabel, ]
  SCdata.diw[ ,SClabel] <- SCdata.sum.merge[ ,SClabel] / plotdata.tmp$fit[1]
}
saveRDS(SCdata.diw, paste0(interfileFolder, "/SCdata.diw_SA", ds.resolution,"CV", CVthr, "_validat_', addcovariate, '.rds"))
if (str_detect(wdpath, "cuizaixu_lab")){
  covariates<-paste0("sex+mean_fd+", addcovariate)
  dataname<-"SCdata.diw"
  smooth_var<-"age"
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-names(SCdata.sum.merge)[1+x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    return(gamresult)
  }, mc.cores = 50)
  
  saveRDS(resultsum, paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE_validat_', addcovariate, '.rds'))
}else{resultsum <- readRDS(paste0(interfileFolder, '/gammodel', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE_validat_', addcovariate, '.rds'))}

# gam results
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-names(SCdata.sum.merge)[1+x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = FALSE, mod_only=FALSE)
    gamresult<-as.data.frame(gamresult)
    return(gamresult)
  }, mc.cores = 50)
  gamresultsum.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
  gamresultsum.df[,c(2:18)]<-lapply(gamresultsum.df[,c(2:18)], as.numeric)
  summary(gamresultsum.df)
  saveRDS(gamresultsum.df, paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE_validat_', addcovariate, '.rds'))
}else{gamresultsum.df <- readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE_validat_', addcovariate, '.rds'))}

# outputs
#1. gamresults', elementnum, '_over8_validat_totalstrength.rds contains statistic variables from gam models
#2. gammodel', elementnum, '_over8_validat_totalstrength.rds contains elementnum gam model files
#3. plotdatasum.df_SA', ds.resolution,'_sumSCinvnode_validat_totalstrength.rds raw plot data
#4. gammodel', elementnum, '_sumSCinvnode_over8_CV75_scale_TRUE_validat_totalstrength.rds the developmental model of the de-weighted SC strength.
# SC strength were de-weighted by the initial SC strength at age of 8 estimated through gam model.
#5. gamresults', elementnum, '_sumSCinvnode_over8_CV75_scale_TRUE_validat_totalstrength.rds the gam stats estimated using scaled data. 
# Mean 2nd derivatives should be calculated using scaled data.

