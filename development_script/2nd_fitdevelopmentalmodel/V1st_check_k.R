# This script is for check the k value.
rm(list=ls())
library(mgcv)
library(parallel)
interfileFolder_HCP <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
interfileFolder_ABCD <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
functionFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
source(paste0(functionFolder, '/gamsmooth.R'))
source(paste0(functionFolder, '/gammsmooth.R'))
SCdata.hcp <- readRDS(paste0(interfileFolder_HCP, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))
SCdata.ABCD <- readRDS(paste0(interfileFolder_ABCD, '/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds'))
# load models
gammodel.HCP.k3 <- readRDS(paste0(interfileFolder_HCP, '/gammodel78_sumSCinvnode_over8_CV75.rds'))
gammodel.ABCD.k3 <- readRDS(paste0(interfileFolder_ABCD, '/gammodel78_sumSCinvnode_over8_siteall_CV75.rds'))

# HCPD k=3~6
AIC_k3_6 <- data.frame(region=paste0("SC.", 1:78), AIC.k3=rep(NA,78), AIC.k4=rep(NA,78),
                       AIC.k5=rep(NA,78),AIC.k6=rep(NA,78))
SCdata.hcp$sex <- as.factor(SCdata.hcp$sex)
covariates<-"sex+mean_fd"
dataname<-"SCdata.hcp"
smooth_var<-"age"
for (k in c(4,5,6)){
  for (i in 1:78){
    model.k3.tmp <- gammodel.HCP.k3[[i]]
    AIC_k3_6$region[i] <- as.character(model.k3.tmp$terms[[2]])
    region <- AIC_k3_6$region[i]
    model.kN.tmp<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=k, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    df.tmp <- AIC(model.k3.tmp, model.kN.tmp)
    AIC_k3_6$AIC.k3[i] <- df.tmp$AIC[1]
    AIC_k3_6[i, paste0('AIC.k', k)] <- df.tmp$AIC[2]
  }
}

min.index <- apply(AIC_k3_6, 1, function(x) which.min(x))
AIC_k3_6$min.index <- min.index+1
table(AIC_k3_6$min.index) # 52/78 edges have the minimal AIC when k=3.

# ABCD k=3~6
AIC_k3_6.abcd <- data.frame(region=paste0("SC.", 1:78), AIC.k3=rep(NA,78), AIC.k4=rep(NA,78),
                       AIC.k5=rep(NA,78),AIC.k6=rep(NA,78))
SCdata.ABCD$sex <- as.factor(SCdata.ABCD$sex)
covariates<-"sex+mean_fd"
dataname<-"SCdata.ABCD"
smooth_var<-"age"
for (k in c(4,5,6)){
  for (i in 1:78){
    model.k3.tmp <- gammodel.ABCD.k3[[i]]
    AIC_k3_6.abcd$region[i] <- as.character(model.k3.tmp$gam$terms[[2]])
    region <- AIC_k3_6.abcd$region[i]
    model.kN.tmp<-gamm.fit.smooth(region, dataname, smooth_var, covariates, knots=k, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    df.tmp <- AIC(model.k3.tmp$mer, model.kN.tmp$mer)
    AIC_k3_6.abcd$AIC.k3[i] <- df.tmp$AIC[1]
    AIC_k3_6.abcd[i, paste0('AIC.k', k)] <- df.tmp$AIC[2]
  }
}

min.index <- apply(AIC_k3_6.abcd, 1, function(x) which.min(x))
AIC_k3_6.abcd$min.index <- min.index+1
table(AIC_k3_6.abcd$min.index) # 71/78 edges have the minimal AIC when k=3.



