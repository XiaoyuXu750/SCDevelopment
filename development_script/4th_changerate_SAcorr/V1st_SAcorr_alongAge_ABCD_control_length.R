## This script is to generate the alignment between developmental effect and connectional axis.
## Derivatives and posterior derivatives generated from the 4th step will be used.
## Derivatives will be used to draw 
## Posterior derivatives will be used to draw 
## Spearman correlation analysis will be used to estimate the correlation coefficients.
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
resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA12'

#### load data
CVthr = 75
derivative.posterior.df <- readRDS(paste0(resultFolder, '/derivative.posterior.df.SA12_CV', CVthr,'.rds'))
derivative.df <- readRDS(paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
mean_length <- readRDS(paste0(interfileFolder, '/mean_length.rds'))
if (CVthr==75){
  mean_length.df <- mean_length[[1]]
  meanlength <- colMeans(mean_length.df[,1:78], na.rm=T)
}else{mean_length.df <- mean_length[[2]]; meanlength <- colMeans(mean_length.df[,1:78], na.rm=T)}
summary(meanlength)
#### source function
source(paste0(functionFolder, '/gamderivatives.R'))
source(paste0(functionFolder, "/SCrankcorr_beforegam.R"))
source(paste0(functionFolder, '/colorbarvalue.R'))
perm.id.full<-readRDS("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA12_sphericalrotations_N10000.rds")

#### calculate SC rank in 12*12
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

#### calculate correlation between posterior derivatives and SC rank
deri.SCrank.posterior.diw.corr <- data.frame(matrix(NA, 1000, 1000))
rownames(deri.SCrank.posterior.diw.corr) <- paste0("draw.", c(1:1000))
colnames(deri.SCrank.posterior.diw.corr) <- paste0("age.", c(1:1000))
# drawtime is the time of draws
# diw is the short of "divide by max SC weight"
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
    deri.pos.tmp_control_length <- residuals(lm(deri.pos.tmp~meanlength))
    corr.pos.tmp <- corr.test(deri.pos.tmp_control_length, SCrank12, method = "spearman")$r
    corr.df$corr.pos.tmp[j]<-corr.pos.tmp
  }
  rownames(corr.df) <- paste0("age.", agerange)
  return(corr.df)
}
# drawtime=1000
deri.SCrank.posterior.corr.sum<-mclapply(1:1000, function(x){
  corr.df.tmp <- compute.SC.corr(x)
  return(corr.df.tmp)
}, mc.cores = 4)

deri.SCrank.posterior.corr<-do.call(rbind, lapply(deri.SCrank.posterior.corr.sum, function(x) t(x$corr.pos.tmp)))
deri.SCrank.posterior.corr<-as.data.frame(deri.SCrank.posterior.corr)
write.csv(deri.SCrank.posterior.corr, paste0(resultFolder, '/deri.SCrank12_CV', CVthr,'.posterior.diw.corr_control_length.csv'), row.names = F)

###### extract age of maximal SC rank correlation: posterior median value + 95% CI
agerange <- unique(derivative.posterior.df[[1]]$age)

# diw
age.max.corr.diw <- lapply(c(1:1000), function(x) agerange[median(which.max(round(deri.SCrank.posterior.corr[x,], 4)))])
age.max.corr.diw <- as.numeric(unlist(age.max.corr.diw))
age.max.corr.diw.median <- median(age.max.corr.diw) #median age #bayes
age.max.corr.diw.CI <- quantile(age.max.corr.diw, probs = c(0.025, 0.975)) #compute the credible interval based on all draws
age.max.corr.diw.lower <- age.max.corr.diw.CI[[1]]
age.max.corr.diw.upper <- age.max.corr.diw.CI[[2]]
age.max.corr.diw.median # 13.62856 (CV75); 21.36055 (CV25)
age.max.corr.diw.CI # 13.45854 13.72571 (CV75); 20.59590 21.77764 (CV25)

age.0.corr.diw <- lapply(c(1:1000), function(x) agerange[median(which.min(abs(round(deri.SCrank.posterior.corr[x,], 4)-0)))])
age.0.corr.diw <- as.numeric(unlist(age.0.corr.diw))
age.0.corr.diw.median <- median(age.0.corr.diw) #median age #bayes
age.0.corr.diw.CI <- quantile(age.0.corr.diw, probs = c(0.025, 0.975)) #compute the credible interval based on all draws
age.0.corr.diw.lower <- age.0.corr.diw.CI[[1]]
age.0.corr.diw.upper <- age.0.corr.diw.CI[[2]]
age.0.corr.diw.median # 13.5557 (CV75) ; 16.07747 (CV25)
age.0.corr.diw.CI # 12.70561 13.72571 (CV75) ; 15.72990 16.42504 (CV25)
#### median corr and 95% CI
# diw
posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
# diw 95%CI
posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))

##### plot alignment with SC rank correlation
#############################################
## diw
df.poscorr.diw <- data.frame(age=agerange, median=posterior.corr.diw.median, up.95CI=posterior.corr.diw.CI$X97.5.,
                             lw.95CI=posterior.corr.diw.CI$X2.5.)
df.poscorr.diw$max.corr.CI <- (df.poscorr.diw$age > age.max.corr.diw.lower & df.poscorr.diw$age < age.max.corr.diw.upper)
df.poscorr.diw$max.corr.window <-df.poscorr.diw$age * df.poscorr.diw$max.corr.CI
df.poscorr.diw$max.corr.window[df.poscorr.diw$max.corr.window==0] <-NA
df.poscorr.diw$zero.corr.CI <- (df.poscorr.diw$age > age.0.corr.diw.lower & df.poscorr.diw$age < age.0.corr.diw.upper)
df.poscorr.diw$zero.corr.window <-df.poscorr.diw$age * df.poscorr.diw$zero.corr.CI
df.poscorr.diw$zero.corr.window[df.poscorr.diw$zero.corr.window==0] <-NA
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
  labs(x="Age (years)", y="Alignment with\nconnectional axis (rho)")+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),
        plot.title = element_text(size=15, hjust = 0.5),
        axis.line = element_line(linewidth = 0.5),
        aspect.ratio =0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15))
ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/SA12_posDeriv_divweight_corr_controllength.tiff'), width=12, height=12, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/SA12_posDeriv_divweight_corr_controllength.svg'), width=14, height=13, units="cm")


