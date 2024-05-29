## Validation for different matrix resolution.
## This script is to generate the alignment between developmental effect and S-A connectional axis.
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
ds.resolution <- 7
elementnum <- ds.resolution*(ds.resolution+1) /2
# set path
wdpath <- getwd()
if (str_detect(wdpath, "Users")){
  resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA', ds.resolution)
}else{
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_ABCD'
}

#### load data
CVthr = 75
derivative.posterior.df <- readRDS(paste0(resultFolder, '/derivative.posterior.df.SA', ds.resolution, '_CV', CVthr,'.rds'))
derivative.df <- readRDS(paste0(resultFolder, '/derivative.df', elementnum, '_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, "/SCrankcorr.R"))
source(paste0(functionFolder, '/colorbarvalue.R'))
#### calculate S-A connectional axis in ds.resolution*ds.resolution
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

#### calculate correlation between posterior derivatives and S-A connectional axis
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
  # estimate rho at 1,000 age points
  for (j in 1:1000){
    deri.pos.tmp <- deriv.SAds.resolution.drawtmp$deri.pos[deriv.SAds.resolution.drawtmp$age==agerange[j]]
    corr.pos.tmp <- corr.test(deri.pos.tmp, SCrankds.resolution, method = "spearman")$r
    corr.df$corr.pos.tmp[j]<-corr.pos.tmp
  }
  rownames(corr.df) <- paste0("age.", agerange)
  return(corr.df)
}

# compute correlation coefficients between S-A connectional axis and 1,000 posterior derivatives at 1,000 age points.
if (str_detect(wdpath, "cuizaixu_lab")){
  deri.SCrank.posterior.corr.sum<-mclapply(1:1000, function(x){
    corr.df.tmp <- compute.SC.corr(x)
    return(corr.df.tmp)
  }, mc.cores = 40)
  
  deri.SCrank.posterior.corr<-do.call(rbind, lapply(deri.SCrank.posterior.corr.sum, function(x) t(x$corr.pos.tmp)))
  deri.SCrank.posterior.corr<-as.data.frame(deri.SCrank.posterior.corr)
  write.csv(deri.SCrank.posterior.corr, paste0(resultFolder, '/deri.SCrank', ds.resolution, '_CV', CVthr,'.posterior.diw.corr.csv'), row.names = F)
}

###### extract age of maximal / zero S-A connectional axis correlation: posterior median value + 95% CI
deri.SCrank.posterior.corr <- read.csv(paste0(resultFolder, '/deri.SCrank', ds.resolution, '_CV', CVthr,'.posterior.diw.corr.csv'))
agerange <- unique(derivative.posterior.df[[1]]$age)
#### median corr and 95% CI
posterior.corr.diw.median <- lapply(c(1:1000), function(x) median(round(deri.SCrank.posterior.corr[,x],4)))
posterior.corr.diw.median <- as.numeric(unlist(posterior.corr.diw.median))
posterior.corr.diw.CI <- lapply(c(1:1000), function(x) quantile(round(deri.SCrank.posterior.corr[,x],4), probs=c(0.025, 0.975)))
posterior.corr.diw.CI <- do.call(rbind, lapply(posterior.corr.diw.CI, function(x) data.frame(t(x))))

##### plot alignment with S-A connectional axis correlation
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
  scale_y_continuous(breaks=c(-0.8, -0.4, 0.0))+
  theme_classic()+
  labs(x="Age (years)", y="Alignment with \nS-A connectional axis (rho)")+
  theme(axis.text=element_text(size=24, color="black"), 
                    axis.title =element_text(size=24),
                    plot.title = element_text(size=15, hjust = 0.5),
                    axis.line = element_line(linewidth = 0.55),
                    axis.ticks = element_line(linewidth = 0.55),
                    aspect.ratio =1.15,
                    plot.background=element_rect(fill="transparent"),
                    panel.background=element_rect(fill="transparent"),plot.margin = margin(t=10, r=0, b=0, l=0),
                    legend.title=element_text(size=15),
                    legend.text=element_text(size=15))
ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/SA', ds.resolution, '_posDeriv_divweight_corr.tiff'), width=12, height=12, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/SA', ds.resolution, '_posDeriv_divweight_corr.svg'), width=13.5, height=15, units="cm")

### correlation plot at different age
#############################################
agerange <- unique(derivative.df$age)
df.poscorr.diw$median[df.poscorr.diw$age==min(agerange)]
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
#SA17: r=-0.6450967 ; SA7: r=-0.8734091

# age=13.7
ntmp<-which.min(abs(agerange-max(agerange))) 
agerange[ntmp]
df.age13 <-as.data.frame(matrix(NA, nrow=1, ncol=9))
names(df.age13)<-names(derivative.df)
for (i in 1:elementnum){
  df.tmp <- derivative.df[derivative.df$label_ID==paste0("SC.", i, "_h"),]
  df.tmp <- df.tmp[df.tmp$age==agerange[ntmp], ]
  df.age13 <- rbind(df.age13, df.tmp)
}
df.age13 <- df.age13[-1,]
SCrankcorr(df.age13, "derivative", ds.resolution, perm.id.full) 
#SA17: r=-0.07378852 ; SA7: r=-0.1831121 

## display 2 lines in one plot
df.agemerge <- rbind(df.age8, df.age13)
df.agemerge$age <- as.factor(df.agemerge$age)
ggplot(data=df.agemerge)+
  geom_smooth(aes(x=rep(SCrankds.resolution, 2), y=derivative, group=age), linewidth=1.2,color="black", method ="lm", 
              se=T)+
  scale_x_continuous(breaks=c(0,40,80))+
  labs(x="S-A connectional axis rank", y="SC change rate")+
  theme_classic()+
  theme(axis.text=element_text(size=21, color="black"), 
        axis.title =element_text(size=21),aspect.ratio = 1,
        plot.title = element_text(size=15, hjust = 0.1, vjust=-5),
        legend.position = "none", plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/deri.diw_corr_SCrank', ds.resolution, 'ageAll.tiff'), width = 14, height = 12, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/deri.diw_corr_SCrank', ds.resolution, 'ageAll.svg'), dpi=600, width = 14, height = 12, units="cm")
#############################################

### derivative plots
##################################
SCrank <- data.frame(label_ID=paste0("SC.",c(1:elementnum), "_h"), SCrankds.resolution=SCrankds.resolution)
SCrank$SCrankds.resolution <- rank(SCrank$SCrankds.resolution, ties.method = "first")
derivative.df.merge <- merge(derivative.df, SCrank, by="label_ID", all.x=T)

## line plot to present change rate
ggplot(data=derivative.df.merge)+
  geom_line(aes(x=age, y=derivative, group=label_ID, color=SCrankds.resolution),size=1.4, alpha=1)+
  scale_color_distiller(type="seq", palette = "RdBu")+
  scale_fill_distiller(type="seq", palette = "RdBu")+
  #geom_vline(xintercept = c(8.92, 13.75), linetype=2, size=0.5)+
  scale_y_continuous(limits=c(-0.076, 0.065), breaks=c(-0.06, -0.03, 0, 0.03, 0.06), labels=c(-6,-3,0,3,6))+
  labs(x="Age (years)", y="SC change rate", color="Axis rank")+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 1.1,
        plot.title = element_text(size=15, hjust = 0.5),
        axis.line = element_line(linewidth=0.5),axis.ticks = element_line(linewidth=0.5),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/derivative_diw_SA', ds.resolution, '_changerate.tiff'), width = 20, height = 16, units="cm")
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/derivative_diw_SA', ds.resolution, '_changerate.svg'), width = 13, height = 13, units="cm")
write.csv(derivative.df.merge, paste0(interfileFolder, '/derivative_df_', ds.resolution, '_merge_CV', CVthr,'.csv'), row.names = F)









