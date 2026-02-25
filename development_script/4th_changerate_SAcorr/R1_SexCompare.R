library(mgcv)
library(parallel)
library(tidyverse)
wdpath <- getwd()
# set resolution
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
# set path
if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_HCPD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-'/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Figure_HCPD_final/SA12'
}else{
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'
}

# set which consistency threshold used to filter spurious streamlines
CVthr=75
# load data
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata_SA", ds.resolution,"_CV75_sumSCinvnode.sum.msmtcsd.combatgam.rds"))
nrow(SCdata.sum.merge)
SCdata.sum.merge$sex <- as.factor(SCdata.sum.merge$sex)
# source function
source(paste0(functionFolder, '/gamsmooth.R'))
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/gamderivatives.R'))
source(paste0(functionFolder, '/SCrankcorr.R'))
detectCores()

SCrank.df <- SCrankcorr(data.frame(rand=rnorm(78)), "rand", 12, T)
SCrank12<-SCrank.df$SCrank

# Define function
compute.SC.corr <- function(drawtime){
  deriv.SA12.drawtmp <- data.frame(age=rep(NA, 78*1000), deri.pos=rep(NA, 78*1000),
                                   SClabel=rep(NA, 78*1000))
  for (i in 1:78){
    df.tmp <- derivative.posterior.df[[i]]
    df.tmp <- df.tmp[df.tmp$draw==paste0("draw", drawtime),]
    lwth <- (i-1)*1000 +1
    upth <- i*1000
    deriv.SA12.drawtmp$age[lwth:upth]<-df.tmp$age
    deriv.SA12.drawtmp$deri.pos[lwth:upth]<-df.tmp$posterior.derivative
    deriv.SA12.drawtmp$SClabel[lwth:upth]<-paste0("SC.", i)
  }
  agerange <- deriv.SA12.drawtmp$age[1:1000]
  corr.df <- data.frame(corr.pos.tmp=rep(NA,1000))
  for (j in 1:1000){
    deri.pos.tmp <- deriv.SA12.drawtmp$deri.pos[deriv.SA12.drawtmp$age==agerange[j]]
    corr.pos.tmp <- corr.test(deri.pos.tmp, SCrank12, method = "spearman")$r
    corr.df$corr.pos.tmp[j]<-corr.pos.tmp
  }
  rownames(corr.df) <- paste0("age.", agerange)
  return(corr.df)
}

DerivativeComp <- function(dataname, ModReturn=F){
  SCdata.tmp <- get(dataname)
  covariates<-"mean_fd"
  smooth_var<-"age"
  # 1) fit models
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-grep("SC.", names(SCdata.sum.merge), value=T)[x]
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=TRUE, stats_only = TRUE, mod_only=TRUE)
    return(gamresult)
  }, mc.cores = 5)
  # 2) estimate plot data
  gammodelsum <- resultsum
  agevect <- SCdata.sum.merge$age
  plotdatasum<-mclapply(1:elementnum, function(x){
    modobj<-gammodelsum[[x]]
    if (is.na(modobj[1])){
      plotdata<-NA
    }else{
      plotdata<- plotdata_generate(modobj, NA, smooth_var="age")
      plotdata$SC_label <- names(plotdata)[13]
      plotdata[,13] <- NULL
    }
    return(plotdata)
  }, mc.cores = 5)
  plotdatasum.df <- do.call(rbind, lapply(plotdatasum, function(x) data.frame(x)))
  # 3) De-weighted data
  SCdata.diw <- SCdata.tmp
  for (i in 1:elementnum){
    SClabel <- paste0("SC.", i, "_h")
    plotdata.tmp <- plotdatasum.df[plotdatasum.df$SC_label==SClabel, ]
    SCdata.diw[ ,SClabel] <- SCdata.tmp[ ,SClabel] / plotdata.tmp$fit[1]
  }
  assign("SCdata.diw", SCdata.diw, envir = .GlobalEnv)
  # 4) fit final model
  dataname<-"SCdata.diw"
  resultsum <- mclapply(1:elementnum, function(x){
    SClabel<-paste0("SC.", x, "_h")
    region<-SClabel
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, 
                              set_fx=T, stats_only = TRUE, mod_only=TRUE)
    return(gamresult)
  }, mc.cores = 5)
  if (ModReturn == T){
    gamresult <- mclapply(1:elementnum, function(x){
      SClabel<-paste0("SC.", x, "_h")
      region<-SClabel
      gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, 
                                set_fx=T, stats_only = TRUE, mod_only=F)
      return(gamresult)
    }, mc.cores = 50)
  }

  # 5) Calculate derivatives
  gammodelsum <- resultsum
  derivative.posterior.df<-mclapply(1:elementnum, function(x){
    SClabel.tmp <- paste0("SC.", x, "_h")
    modobj<-gammodelsum[[x]]
    draws<-1000
    increments<-1000
    derivdata<- gam.derivatives(modobj, "age",smoothvector=NA, draws, increments, 
                                return_posterior_derivatives = TRUE)
    derivdata$SCrank<-SCrank12[x]
    meanSC <- mean(modobj$model[,SClabel.tmp],na.rm=T)
    derivdata$meanSC<-meanSC
    maxSC <- max(modobj$model[,SClabel.tmp],na.rm=T)
    derivdata$maxSC<-maxSC
    return(derivdata)
  }, mc.cores = 5)
 
  if (ModReturn == F){
    return(derivative.posterior.df)
  }else{
    return(list(gammodelsum, derivative.posterior.df, gamresult))
  }
  
}

ComputeAge0 <- function(derivative.posterior.df, PlotdfReturn=F){
  # 6) Estimate age window
  deri.SCrank.posterior.diw.corr <- data.frame(matrix(NA, 1000, 1000))
  rownames(deri.SCrank.posterior.diw.corr) <- paste0("draw.", c(1:1000))
  colnames(deri.SCrank.posterior.diw.corr) <- paste0("age.", c(1:1000))
  
  deri.SCrank.posterior.corr.sum<-mclapply(1:1000, function(x){
    corr.df.tmp <- compute.SC.corr(x)
    return(corr.df.tmp)
  }, mc.cores = 5)
  deri.SCrank.posterior.corr<-do.call(rbind, lapply(deri.SCrank.posterior.corr.sum, function(x) t(x$corr.pos.tmp)))
  deri.SCrank.posterior.corr<-as.data.frame(deri.SCrank.posterior.corr)
  
  # 7) Compute across zero age
  agerange <- unique(derivative.posterior.df[[1]]$age)
  
  # 1. age window when alignment is equal to zero
  age.0.corr.diw <- lapply(c(1:1000), function(x) agerange[median(which.min(abs(round(deri.SCrank.posterior.corr[x,], 4)-0)))])
  age.0.corr.diw <- as.numeric(unlist(age.0.corr.diw))
  age.0.corr.diw.median <- median(age.0.corr.diw) #median age #bayes
  age.0.corr.diw.CI <- quantile(age.0.corr.diw, probs = c(0.025, 0.975))
  age.0.corr.diw.lower <- age.0.corr.diw.CI[[1]]
  age.0.corr.diw.upper <- age.0.corr.diw.CI[[2]]
  
  # 8) Plot data
  ##### plot alignment with S-A connectional axis correlation (Figure 3 (B))
  # diw
  posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
  posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
  # diw 95%CI
  posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
  posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))
  df.poscorr.diw <- data.frame(age=agerange, median=posterior.corr.diw.median, up.95CI=posterior.corr.diw.CI$X97.5.,
                               lw.95CI=posterior.corr.diw.CI$X2.5.)
  df.poscorr.diw$zero.corr.CI <- (df.poscorr.diw$age > age.0.corr.diw.lower & df.poscorr.diw$age < age.0.corr.diw.upper)
  df.poscorr.diw$zero.corr.window <-df.poscorr.diw$age * df.poscorr.diw$zero.corr.CI
  df.poscorr.diw$zero.corr.window[df.poscorr.diw$zero.corr.window==0] <-NA
  loess.median <- loess(median~age, data=df.poscorr.diw, span=0.2)
  loess.lw <- loess(lw.95CI~age, data=df.poscorr.diw, span=0.2)
  loess.up <- loess(up.95CI~age, data=df.poscorr.diw, span=0.2)
  df.poscorr.diw$median.loess <- loess.median$fitted
  df.poscorr.diw$lw.95CI.loess <- loess.lw$fitted
  df.poscorr.diw$up.95CI.loess <- loess.up$fitted
  
  age.0 <- c(age.0.corr.diw.median, age.0.corr.diw.lower, age.0.corr.diw.upper)
  if (PlotdfReturn==F){
    return(age.0)
  }else{
    return(list(age.0, df.poscorr.diw, age.0.corr.diw))
  }
}


n <- commandArgs(trailingOnly = TRUE)
n <- as.numeric(n)

if (n == 0){
  SCdata.F <- SCdata.sum.merge %>% filter(sex =="F")
  SCdata.M <- SCdata.sum.merge %>% filter(sex =="M")

  derivative.posterior.df <- DerivativeComp("SCdata.F")
  age.0.F <- ComputeAge0(derivative.posterior.df)

  derivative.posterior.df <- DerivativeComp("SCdata.M")
  age.0.M <- ComputeAge0(derivative.posterior.df)

}else{
  set.seed(925+n)
  SCdata.rand <- SCdata.sum.merge %>%
    mutate(sex2 = sex[sample(1:nrow(SCdata.sum.merge),nrow(SCdata.sum.merge), F)])
  SCdata.F <- SCdata.rand %>% filter(sex2 =="F")
  SCdata.M <- SCdata.rand %>% filter(sex2 =="M")

  derivative.posterior.df <- DerivativeComp("SCdata.F")
  age.0.F <- ComputeAge0(derivative.posterior.df)

  derivative.posterior.df <- DerivativeComp("SCdata.M")
  age.0.M <- ComputeAge0(derivative.posterior.df)
}

if (! dir.exists(paste0(interfileFolder, "/SexCompare_Permut"))){
  dir.create(paste0(interfileFolder, "/SexCompare_Permut"))
}

df <- cbind(t(age.0.F), t(age.0.M))
df <- as.data.frame(df)
names(df) <- c("Median.F", "CI.low.F", "CI.up.F", "Median.M", "CI.low.M", "CI.up.M")

write.csv(df, paste0(interfileFolder, "/SexCompare_Permut/ZeroAlign_By_Sex_", n, ".csv"), row.names = F)
