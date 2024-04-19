## Validation for different matrix resolution
## This script is to generate the alignment between developmental effect and connectional axis.
## Derivatives and posterior derivatives generated from the 2nd step will be used.
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
ds.resolution <- 17
elementnum <- ds.resolution*(ds.resolution+1) /2
# set path
wdpath <- getwd()
if (str_detect(wdpath, "Users")){
  resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_HCPD'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final/SA', ds.resolution)
}else{
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_HCPD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_HCPD'
}

#### load data
CVthr = 75
derivative.posterior.df <- readRDS(paste0(resultFolder, '/derivative.posterior.df.SA', ds.resolution, '_CV', CVthr,'.rds'))
derivative.df <- readRDS(paste0(resultFolder, '/derivative.df', elementnum, '_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/gamderivatives.R'))
source(paste0(functionFolder, "/SCrankcorr.R"))
source(paste0(functionFolder, '/colorbarvalue.R'))
perm.id.full<-readRDS(paste0("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA", ds.resolution, "_sphericalrotations_N10000.rds"))

#### calculate SC rank in ds.resolution*ds.resolution
########################
Matrixds.resolution<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
indexupds.resolution <- upper.tri(Matrixds.resolution)
indexsaveds.resolution <- !indexupds.resolution
Matrixds.resolution.index<-Matrixds.resolution
#(ds.resolution+1)*ds.resolution/2=elementnum
Matrixds.resolution.index[indexsaveds.resolution]<-c(1:elementnum)
Matrixds.resolution.SCrank<-Matrixds.resolution
for (x in 1:ds.resolution){
  for (y in 1:ds.resolution){
    Matrixds.resolution.SCrank[x,y]<-x^2+y^2
  }
}
Matrixds.resolution.SCrank[indexupds.resolution]<-NA
SCrankds.resolution<-rank(Matrixds.resolution.SCrank[indexsaveds.resolution], ties.method = "average")

#### calculate correlation between posterior derivatives and SC rank
deri.SCrank.posterior.diw.corr <- data.frame(matrix(NA, 1000, 1000))
rownames(deri.SCrank.posterior.diw.corr) <- paste0("draw.", c(1:1000))
colnames(deri.SCrank.posterior.diw.corr) <- paste0("age.", c(1:1000))
# drawtime is the time of draws
compute.SC.corr <- function(drawtime){
  deriv.SAds.resolution.drawtmp <- data.frame(age=rep(NA, elementnum*1000), deri.pos=rep(NA, elementnum*1000),
                                   SClabel=rep(NA, elementnum*1000))
  for (i in 1:elementnum){
    df.tmp <- derivative.posterior.df[[i]]
    df.tmp <- df.tmp[df.tmp$draw==paste0("draw", drawtime),]
    lwth <- (i-1)*1000 +1
    upth <- i*1000
    deriv.SAds.resolution.drawtmp$age[lwth:upth]<-df.tmp$age
    deriv.SAds.resolution.drawtmp$deri.pos[lwth:upth]<-df.tmp$posterior.derivative
    deriv.SAds.resolution.drawtmp$SClabel[lwth:upth]<-paste0("SC.", i)
  }
  agerange <- deriv.SAds.resolution.drawtmp$age[1:1000]
  corr.df <- data.frame(corr.pos.tmp=rep(NA,1000))
  for (j in 1:1000){
    deri.pos.tmp <- deriv.SAds.resolution.drawtmp$deri.pos[deriv.SAds.resolution.drawtmp$age==agerange[j]]
    corr.pos.tmp <- corr.test(deri.pos.tmp, SCrankds.resolution, method = "spearman")$r
    corr.df$corr.pos.tmp[j]<-corr.pos.tmp
  }
  rownames(corr.df) <- paste0("age.", agerange)
  return(corr.df)
}

# compute correlation coefficients between connectional axis and 1,000 posterior derivatives at 1,000 age points.
if (str_detect(wdpath, "cuizaixu_lab")){
  deri.SCrank.posterior.corr.sum<-mclapply(1:1000, function(x){
    corr.df.tmp <- compute.SC.corr(x)
    return(corr.df.tmp)
  }, mc.cores = 40)
  
  deri.SCrank.posterior.corr<-do.call(rbind, lapply(deri.SCrank.posterior.corr.sum, function(x) t(x$corr.pos.tmp)))
  deri.SCrank.posterior.corr<-as.data.frame(deri.SCrank.posterior.corr)
  write.csv(deri.SCrank.posterior.corr, paste0(resultFolder, '/deri.SCrank', ds.resolution, '_CV', CVthr,'.posterior.diw.corr.csv'), row.names = F)
}

###### extract age of maximal / zero SC rank correlation: posterior median value + 95% CI
deri.SCrank.posterior.corr <- read.csv(paste0(resultFolder, '/deri.SCrank', ds.resolution, '_CV', CVthr,'.posterior.diw.corr.csv'))
agerange <- unique(derivative.posterior.df[[1]]$age)

# diw
age.max.corr.diw <- lapply(c(1:1000), function(x) agerange[median(which.max(round(deri.SCrank.posterior.corr[x,], 4)))])
age.max.corr.diw <- as.numeric(unlist(age.max.corr.diw))
age.max.corr.diw.median <- median(age.max.corr.diw) #median age #bayes
age.max.corr.diw.CI <- quantile(age.max.corr.diw, probs = c(0.025, 0.975)) #compute the credible interval based on all draws
age.max.corr.diw.lower <- age.max.corr.diw.CI[[1]]
age.max.corr.diw.upper <- age.max.corr.diw.CI[[2]]
age.max.corr.diw.median # 21.4874 (resolution=17); 20.35194 (resolution=7)
age.max.corr.diw.CI # 20.89163 21.80589 (resolution=17); 18.68993 21.47425 (resolution=7)

age.0.corr.diw <- lapply(c(1:1000), function(x) agerange[median(which.min(abs(round(deri.SCrank.posterior.corr[x,], 4)-0)))])
age.0.corr.diw <- as.numeric(unlist(age.0.corr.diw))
age.0.corr.diw.median <- median(age.0.corr.diw) #median age 
age.0.corr.diw.CI <- quantile(age.0.corr.diw, probs = c(0.025, 0.975)) #compute the credible interval based on all draws
age.0.corr.diw.lower <- age.0.corr.diw.CI[[1]]
age.0.corr.diw.upper <- age.0.corr.diw.CI[[2]]
age.0.corr.diw.median # 15.50542 (resolution=17) ; 15.56773 (resolution=7)
age.0.corr.diw.CI # 15.29771 15.75467 (resolution=17) ; 15.15924 15.99007 (resolution=7)
#### median corr and 95% CI
# diw
posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
# diw 95%CI
posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))

##### plot alignment with SC rank correlation (Supplementary Figure 4)
#############################################
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
  geom_line(aes(x=age, y=median.loess), linewidth=1.5)+
  #geom_ribbon(aes(x=max.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#595959")+
  geom_ribbon(aes(x=zero.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#F8B01B", alpha=1)+
  geom_ribbon(aes(x=zero.corr.window, ymin=median.loess-0.045, ymax=median.loess+0.045), fill="#B2182B", alpha=1)+
  #geom_vline(aes(xintercept=age.0.corr.diw.median), color="black",linetype="dashed")+
  theme_classic()+
  labs(x="Age (years)", y="Alignment with connectional axis (rho)")+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22),
        axis.line = element_line(linewidth = 0.5),
        aspect.ratio =0.93,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15))
ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/SA', ds.resolution, '_posDeriv_divweight_corr.tiff'), width=12, height=12, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/SA', ds.resolution, '_posDeriv_divweight_corr.svg'), dpi=600, width=15, height =15, units="cm")

### plot age distribution with 0 corr
ggplot() +
  geom_histogram(aes(age.0.corr.diw, y = ..count..),binwidth = 0.5, linewidth=1,
                 color = "black", fill = "white") +
  geom_vline(aes(xintercept = age.0.corr.diw.median), colour = "red", linetype="solid", linewidth=1)+
  labs(x = NULL, y = NULL) +theme_classic()+
  scale_y_discrete(breaks=NULL)+scale_x_continuous(breaks=NULL)+
  theme(axis.line = element_blank(),
        aspect.ratio = 0.5,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.title = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", size = 20))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Agedistribution_0corr.tiff'), width=6, height = 6, units="cm")
ggsave(paste0(FigureFolder, '/CV', CVthr, '/Alignment_development/Agedistribution_0corr.svg'), width=6, height =6, units="cm")

ggplot() +
  geom_histogram(aes(age.max.corr.diw, y = ..count..),binwidth = 0.5, linewidth=0.8,
                 color = "black", fill = "white") +
  geom_vline(aes(xintercept = age.max.corr.diw.median), colour = "red", linetype="solid")+
  labs(x = NULL, y = NULL) +theme_classic()+
  scale_y_discrete(breaks=NULL)+scale_x_continuous(breaks=NULL)+
  theme(axis.line = element_blank(),
        aspect.ratio = 0.5,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.title = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", size = 20))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Agedistribution_maxcorr.tiff'), width=6, height = 6, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Agedistribution_maxcorr.svg'), width=6, height =6, units="cm")

### correlation plot at different age (Figure 4 (d))
#############################################
agerange <- unique(derivative.df$age)
df.poscorr.diw$median[df.poscorr.diw$age==min(agerange)] #-0.5738
# age = 8.083333
min(agerange)
df.age8 <-as.data.frame(matrix(NA, nrow=1, ncol=9))
names(df.age8)<-names(derivative.df)
for (i in 1:elementnum){
  df.tmp <- derivative.df[derivative.df$label_ID==paste0("SC.", i, "_h"),]
  df.tmp <- df.tmp[df.tmp$age==min(agerange), ]
  df.age8 <- rbind(df.age8, df.tmp)
}
df.age8 <- df.age8[-1,]

SCrankcorr(df.age8, "derivative", ds.resolution, perm.id.full) 
#SA17: r=-0.6387991  p=0.010825; SA7: r=-0.7850007  p=0.02305

# age = age.0.corr.diw.median
ntmp<-which.min(abs(agerange-age.0.corr.diw.median))
agerange[ntmp]
df.poscorr.diw$median[round(df.poscorr.diw$age,2)==round(age.0.corr.diw.median,2)] #0.013
df.age15 <-as.data.frame(matrix(NA, nrow=1, ncol=9))
names(df.age15)<-names(derivative.df)
for (i in 1:elementnum){
  df.tmp <- derivative.df[derivative.df$label_ID==paste0("SC.", i, "_h"),]
  df.tmp <- df.tmp[df.tmp$age==agerange[ntmp], ]
  df.age15 <- rbind(df.age15, df.tmp)
}
df.age15 <- df.age15[-1,]
SCrankcorr(df.age15, "derivative", ds.resolution, perm.id.full) 
#SA17: r=-0.003032359 p=0.4393; SA7: r=-0.02381278 p=0.390325

# age = 21.92
ntmp<-which.min(abs(agerange-max(agerange)))
agerange[ntmp]
df.poscorr.diw$median[round(df.poscorr.diw$age,2)==round(age.max.corr.diw.median)] #0.641
df.age21 <-as.data.frame(matrix(NA, nrow=1, ncol=9))
names(df.age21)<-names(derivative.df)
for (i in 1:elementnum){
  df.tmp <- derivative.df[derivative.df$label_ID==paste0("SC.", i, "_h"),]
  df.tmp <- df.tmp[df.tmp$age==agerange[ntmp], ]
  df.age21 <- rbind(df.age21, df.tmp)
}
df.age21 <- df.age21[-1,]
SCrankcorr(df.age21, "derivative", ds.resolution, perm.id.full) 
#SA17: r=0.7187529 p=0.00025; SA7: r=0.8411113 p=0.0059

## display 3 lines in one plot
df.agemerge <- rbind(df.age8, df.age15, df.age21)
df.agemerge$age <- as.factor(df.agemerge$age)
ggplot(data=df.agemerge)+
  geom_smooth(aes(x=rep(SCrankds.resolution, 3), y=derivative, group=age),color="black", method ="lm", 
              se=T)+
  labs(x="Connectional axis rank", y="SC change rate")+
  theme_classic()+
  theme(axis.text=element_text(size=15.4, color="black"), 
        axis.title =element_text(size=15.4),aspect.ratio = 0.65,
        plot.title = element_text(size=15, hjust = 0.1, vjust=-5),axis.line = element_line(linewidth = 0.35),
        axis.ticks  = element_line(linewidth = 0.35),
        legend.position = "none", plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/deri.diw_corr_SCrank', ds.resolution, 'ageAll.tiff'), width = 14, height = 12, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/deri.diw_corr_SCrank', ds.resolution, 'ageAll.svg'), width = 12, height = 10, units="cm")
#############################################

### derivative plots (Supplementary Figure4)
##################################
SCrank <- data.frame(label_ID=paste0("SC.",c(1:elementnum), "_h"), SCrankds.resolution=SCrankds.resolution)
SCrank$SCrankds.resolution <- rank(SCrank$SCrankds.resolution, ties.method = "first")
derivative.df.merge <- merge(derivative.df, SCrank, by="label_ID", all.x=T)

## line plot to present change rate
ggplot(data=derivative.df.merge)+
  geom_line(aes(x=age, y=derivative, group=label_ID, color=SCrankds.resolution),size=1.4, alpha=1)+
  scale_color_distiller(type="seq", palette = "RdBu")+
  scale_fill_distiller(type="seq", palette = "RdBu")+
  #geom_vline(xintercept = c(8.09, round(age.0.corr.diw.median,1), 21.2), linetype=2, size=0.5)+
  scale_y_continuous(breaks = c(-0.02, 0, 0.02), labels = c("-2", "0", "2"))+
  labs(x="Age (years)", y="SC change rate", color="Axis rank")+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 1,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/derivative_diw_SA', ds.resolution, '_changerate.tiff'), width = 20, height = 16, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/derivative_diw_SA', ds.resolution, '_changerate.svg'), dpi=600, width=14, height =14, units="cm")
write.csv(derivative.df.merge, paste0(interfileFolder, '/derivative_df_', ds.resolution, '_merge_CV', CVthr,'.csv'), row.names = F)
