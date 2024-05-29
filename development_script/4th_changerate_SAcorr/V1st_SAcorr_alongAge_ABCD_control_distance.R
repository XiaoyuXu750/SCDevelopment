## Validation for controling Euclidean distance in ABCD dataset.
## Derivatives and posterior derivatives will be used.

library(tidyverse)
library(R.matlab)
library(psych)
library(gratia)
library(mgcv)
library(parallel)
library(ggplot2)
library(RColorBrewer)
library(reshape)

rm(list = ls())
wdpath <- getwd()
CVthr = 75
if (str_detect(wdpath, "Users")){
  resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA12'
  source(paste0(functionFolder, "/SCrankcorr.R"))

}else if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_ABCD'
}

#### load data
derivative.posterior.df <- readRDS(paste0(resultFolder, '/derivative.posterior.df.SA12_CV', CVthr,'.rds'))
derivative.df <- readRDS(paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
meandistance <- read.csv(paste0(interfileFolder, "/average_EuclideanDistance_12.csv"))
meandistance <- meandistance$Edistance
#### source function
source(paste0(functionFolder, '/gamderivatives.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))

#### calculate S-A connectional axis in 12*12
########################
Matrix12<-matrix(NA, nrow=12, ncol=12)
indexup12 <- upper.tri(Matrix12)
indexsave12 <- !indexup12
Matrix12.index<-Matrix12
#13*12/2=78
Matrix12.index[indexsave12]<-c(1:78)
#tilt rank
Matrix12.tiltrank<-Matrix12
for (x in 1:12){
  for (y in 1:12){
    Matrix12.tiltrank[x,y]<-(x+y)^2+(x-y)^2
  }
}
Matrix12.tiltrank[indexup12]<-NA
SCrank12<-rank(Matrix12.tiltrank[indexsave12], ties.method = "average")

#### calculate correlation between posterior derivatives and S-A connectional axis
deri.SCrank.posterior.diw.corr <- data.frame(matrix(NA, 1000, 1000))
rownames(deri.SCrank.posterior.diw.corr) <- paste0("draw.", c(1:1000))
colnames(deri.SCrank.posterior.diw.corr) <- paste0("age.", c(1:1000))
# drawtime is the time of draws
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
    deri.pos.tmp_control_distance <- residuals(lm(deri.pos.tmp~meandistance))
    corr.pos.tmp <- corr.test(deri.pos.tmp_control_distance, SCrank12, method = "spearman")$r
    corr.df$corr.pos.tmp[j]<-corr.pos.tmp
  }
  rownames(corr.df) <- paste0("age.", agerange)
  return(corr.df)
}
# drawtime=1000
if (str_detect(wdpath, "cuizaixu_lab")){
  deri.SCrank.posterior.corr.sum<-mclapply(1:1000, function(x){
    corr.df.tmp <- compute.SC.corr(x)
    return(corr.df.tmp)
  }, mc.cores = 50)
  
  deri.SCrank.posterior.corr<-do.call(rbind, lapply(deri.SCrank.posterior.corr.sum, function(x) t(x$corr.pos.tmp)))
  deri.SCrank.posterior.corr<-as.data.frame(deri.SCrank.posterior.corr)
  write.csv(deri.SCrank.posterior.corr, paste0(resultFolder, '/deri.SCrank12_CV', CVthr,'.posterior.diw.corr_control_distance.csv'), row.names = F)
  
}else{
  deri.SCrank.posterior.corr <- read.csv(paste0(resultFolder, '/deri.SCrank12_CV', CVthr,'.posterior.diw.corr_control_distance.csv'))
}

###### extract age of maximal S-A connectional axis: posterior median value + 95% CI
agerange <- unique(derivative.posterior.df[[1]]$age)
posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))

##### plot alignment with S-A connectional axis
#############################################
df.poscorr.diw <- data.frame(age=agerange, median=posterior.corr.diw.median, up.95CI=posterior.corr.diw.CI$X97.5.,
                             lw.95CI=posterior.corr.diw.CI$X2.5.)
loess.median <- loess(median~age, data=df.poscorr.diw, span=0.2)
loess.lw <- loess(lw.95CI~age, data=df.poscorr.diw, span=0.2)
loess.up <- loess(up.95CI~age, data=df.poscorr.diw, span=0.2)
df.poscorr.diw$median.loess <- loess.median$fitted
df.poscorr.diw$lw.95CI.loess <- loess.lw$fitted
df.poscorr.diw$up.95CI.loess <- loess.up$fitted

ggplot(data=df.poscorr.diw)+
  geom_ribbon(aes(x=age, ymin=lw.95CI.loess, ymax=up.95CI.loess), alpha=0.3)+
  geom_line(aes(x=age, y=median.loess), size=1)+
  #geom_ribbon(aes(x=max.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#595959")+
  #geom_ribbon(aes(x=zero.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#F8B01B", alpha=1)+
  #geom_ribbon(aes(x=zero.corr.window, ymin=median.loess-0.005, ymax=median.loess+0.005), fill="#B2182B", alpha=1)+
  #geom_vline(aes(xintercept=age.0.corr.diw.median), color="black",linetype="dashed")+
  theme_classic()+
  labs(x="Age (years)", y="Alignment with\nS-A connectional axis (rho)")+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),
        plot.title = element_text(size=15, hjust = 0.5),
        axis.line = element_line(linewidth = 0.5),
        aspect.ratio =0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15))
ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/SA12_posDeriv_divweight_corr_controldistance.tiff'), width=12, height=12, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/SA12_posDeriv_divweight_corr_controldistance.svg'), width=14, height=13, units="cm")


