rm(list=ls())
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

## Actual
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
Pjpath = str_split_i(wdpath, "Rcode_SCdevelopment", 1)
source(paste0(Pjpath, "Rcode_SCdevelopment/development_script/2nd_fitdevelopmentalmodel/R4_SexCompare.R"))
detectCores()

SCrank.df <- SCrankcorr(data.frame(rand=rnorm(78)), "rand", 12, T)
SCrank12<-SCrank.df$SCrank


# Actual Female & Male
SCdata.F <- SCdata.sum.merge %>% filter(sex =="F")
SCdata.M <- SCdata.sum.merge %>% filter(sex =="M")
## F
set.seed(925)
Modtmp.F <- DerivativeComp("SCdata.F", T)
gammod.F <- Modtmp.F[[1]]
derivative.posterior.df <- Modtmp.F[[2]]
gamresult.F <- Modtmp.F[[3]]
gamresult.F2 <- do.call(rbind, lapply(gamresult.F, as.data.frame))
gamresult.F2[,2:10] <- lapply(gamresult.F2[,2:10], as.numeric)
summary(gamresult.F2)
saveRDS(gamresult.F2, paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_Female.rds'))

Align0.F <- ComputeAge0(derivative.posterior.df, PlotdfReturn=T)
age.0.F <- Align0.F[[1]]
alignPlotdf.F <- Align0.F[[2]]
age.0.corr.diw.F <- Align0.F[[3]]
saveRDS(alignPlotdf.F, paste0(interfileFolder, '/Alignment_Plotdata_SA12_Female.rds'))
saveRDS(age.0.corr.diw.F, paste0(interfileFolder, '/Alignment_age.0.corr.diw_SA12_Female.rds'))

## M
set.seed(925)
Modtmp.M <- DerivativeComp("SCdata.M", T)
gammod.M <- Modtmp.M[[1]]
derivative.posterior.df <- Modtmp.M[[2]]
gamresult.M <- Modtmp.M[[3]]
gamresult.M2 <- do.call(rbind, lapply(gamresult.M, as.data.frame))
gamresult.M2[,2:10] <- lapply(gamresult.M2[,2:10], as.numeric)
summary(gamresult.M2)
saveRDS(gamresult.M2, paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_Male.rds'))

Align0.M <- ComputeAge0(derivative.posterior.df, PlotdfReturn=T)
age.0.M <- Align0.M[[1]]
alignPlotdf.M <- Align0.M[[2]]
age.0.corr.diw.M <- Align0.M[[3]]
saveRDS(alignPlotdf.M, paste0(interfileFolder, '/Alignment_Plotdata_SA12_Male.rds'))
saveRDS(age.0.corr.diw.M, paste0(interfileFolder, '/Alignment_age.0.corr.diw_SA12_Male.rds'))

## Save out age at zero alignment
df <- cbind(t(age.0.F), t(age.0.M))
df <- as.data.frame(df)
names(df) <- c("Median.F", "CI.low.F", "CI.up.F", "Median.M", "CI.low.M", "CI.up.M")

write.csv(df, paste0(interfileFolder, "/SexCompare_Permut/ZeroAlign_By_Sex_0.csv"), row.names = F)

## 1) Generate plot data
### Female
plotdatasum.F <- mclapply(1:elementnum, function(x){
  modobj<-gammod.F[[x]]
  if (is.na(modobj[1])){
    plotdata<-NA
  }else{
    plotdata<- plotdata_generate(modobj, NA, smooth_var="age")
    plotdata$SC_label <- names(plotdata)[13]
    plotdata[,13] <- NULL
  }
  return(plotdata)
}, mc.cores = 5)
plotdatasum.df.F <- do.call(rbind, lapply(plotdatasum.F, function(x) data.frame(x)))

saveRDS(plotdatasum.df.F, paste0(interfileFolder, "/Fit_plotdf_Female.rds"))

### Male
plotdatasum.M <- mclapply(1:elementnum, function(x){
  modobj<-gammod.M[[x]]
  if (is.na(modobj[1])){
    plotdata<-NA
  }else{
    plotdata<- plotdata_generate(modobj, NA, smooth_var="age")
    plotdata$SC_label <- names(plotdata)[13]
    plotdata[,13] <- NULL
  }
  return(plotdata)
}, mc.cores = 5)
plotdatasum.df.M <- do.call(rbind, lapply(plotdatasum.M, function(x) data.frame(x)))

saveRDS(plotdatasum.df.M, paste0(interfileFolder, "/Fit_plotdf_Male.rds"))

### Draw picture
## Trajectories
Actual.age0 <- read.csv(paste0(interfileFolder, "/SexCompare_Permut/ZeroAlign_By_Sex_0.csv"))
# F
###################
plotdatasum.df.F <- readRDS(paste0(interfileFolder, "/Fit_plotdf_Female.rds"))
gamresult.F <- readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_Female.rds'))

plotdatasum.df.F <- plotdatasum.df.F %>% left_join(gamresult.F, join_by(SC_label == parcel))
lmthr <- max(abs(gamresult.F$partialRsq)) # 0.144
ggplot()+
  geom_line(data=plotdatasum.df.F, aes(x=age, y=fit.ratio, group=SC_label, color=partialRsq), size=1.5, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1, limit=c(-lmthr,lmthr))+
  labs(x="Age (years)", y="SC strength (ratio)")+
  #scale_color_manual(values = rev(brewer.pal(10, "RdBu")))+
  #scale_y_continuous(breaks=c(1.0, 1.1, 1.2))+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=22, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV75/SA12_sumSCinvnode_fit/devcurve_Rsq_fit.ratio_Female.tiff'), width=15, height =14, units = "cm")

alignPlotdf.F <- readRDS(paste0(interfileFolder, '/Alignment_Plotdata_SA12_Female.rds'))
ggplot(data=alignPlotdf.F)+
  geom_ribbon(aes(x=age, ymin=lw.95CI.loess, ymax=up.95CI.loess), alpha=0.3)+
  geom_line(aes(x=age, y=median.loess), linewidth=1.5)+
  #geom_ribbon(aes(x=max.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#595959")+
  geom_ribbon(aes(x=zero.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#F8B01B", alpha=1)+
  geom_ribbon(aes(x=zero.corr.window, ymin=median.loess-0.045, ymax=median.loess+0.045), fill="#B2182B", alpha=1)+
  #geom_vline(aes(xintercept=age.0.corr.diw.median), color="black",linetype="dashed")+
  theme_classic()+
  labs(x="Age (years)", y="rho")+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 0.93,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=22, hjust = 0.5), legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/SA12_posDeriv_divweight_corr_Female.tiff'), width=15, height=14, units="cm")

### plot age distribution with 0 corr
age.0.corr.diw.F <- readRDS(paste0(interfileFolder, '/Alignment_age.0.corr.diw_SA12_Female.rds'))
ggplot() +
  geom_histogram(aes(age.0.corr.diw.F, y = ..count..),binwidth = 0.5, linewidth=1,
                 color = "black", fill = "white") +
  geom_vline(aes(xintercept = Actual.age0$Median.F), colour = "red", linetype="solid", linewidth=1)+
  labs(x = NULL, y = NULL) +theme_classic()+
  scale_y_discrete(breaks=NULL)+scale_x_continuous(breaks=NULL)+
  theme(axis.line = element_blank(),
        aspect.ratio = 0.5,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.title = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", size = 20))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Agedistribution_0corr_Female.tiff'), width=6, height = 6, units="cm")
##################

# M
#################
plotdatasum.df.M <- readRDS(paste0(interfileFolder, "/Fit_plotdf_Male.rds"))
gamresult.M <- readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_Male.rds'))

plotdatasum.df.M <- plotdatasum.df.M %>% left_join(gamresult.M, join_by(SC_label == parcel))
lmthr <- max(abs(gamresult.M$partialRsq)) # 0.129
ggplot()+
  geom_line(data=plotdatasum.df.M, aes(x=age, y=fit.ratio, group=SC_label, color=partialRsq), size=1.5, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1, limit=c(-lmthr,lmthr))+
  labs(x="Age (years)", y="SC strength (ratio)")+
  #scale_color_manual(values = rev(brewer.pal(10, "RdBu")))+
  #scale_y_continuous(breaks=c(1.0, 1.1, 1.2))+
  theme_classic()+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=22, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV75/SA12_sumSCinvnode_fit/devcurve_Rsq_fit.ratio_Male.tiff'), width=15, height =14, units = "cm")

alignPlotdf.M <- readRDS(paste0(interfileFolder, '/Alignment_Plotdata_SA12_Male.rds'))
ggplot(data=alignPlotdf.M)+
  geom_ribbon(aes(x=age, ymin=lw.95CI.loess, ymax=up.95CI.loess), alpha=0.3)+
  geom_line(aes(x=age, y=median.loess), linewidth=1.5)+
  #geom_ribbon(aes(x=max.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#595959")+
  geom_ribbon(aes(x=zero.corr.window, ymin=lw.95CI.loess, ymax=up.95CI.loess), fill="#F8B01B", alpha=1)+
  geom_ribbon(aes(x=zero.corr.window, ymin=median.loess-0.045, ymax=median.loess+0.045), fill="#B2182B", alpha=1)+
  scale_y_continuous(breaks = c(-0.4, 0, 0.4))+
  theme_classic()+
  labs(x="Age (years)", y="rho")+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 0.93,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=22, hjust = 0.5), legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/Alignment_development/SA12_posDeriv_divweight_corr_Male.tiff'), width=15, height=14, units="cm")

### plot age distribution with 0 corr
age.0.corr.diw.M <- readRDS(paste0(interfileFolder, '/Alignment_age.0.corr.diw_SA12_Male.rds'))
ggplot() +
  geom_histogram(aes(age.0.corr.diw.M, y = ..count..),binwidth = 0.5, linewidth=1,
                 color = "black", fill = "white") +
  geom_vline(aes(xintercept = Actual.age0$Median.M), colour = "red", linetype="solid", linewidth=1)+
  labs(x = NULL, y = NULL) +theme_classic()+
  scale_y_discrete(breaks=NULL)+scale_x_continuous(breaks=NULL)+
  theme(axis.line = element_blank(),
        aspect.ratio = 0.5,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        axis.title = element_text(color = "black", size = 15),
        axis.text.x = element_text(color = "black", size = 20))
ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/Agedistribution_0corr_Male.tiff'), width=6, height = 6, units="cm")
#####################

#### Permutation
Actual.age0 <- read.csv(paste0(interfileFolder, "/SexCompare_Permut/ZeroAlign_By_Sex_0.csv"))
# Median.F CI.low.F  CI.up.F Median.M CI.low.M  CI.up.M
# 1  15.1177 14.84076 15.45003 16.17009 15.74082 16.73851

Age0.perm <- list()
for (i in 1:1000){
  tmp <- read.csv(paste0(interfileFolder, "/SexCompare_Permut/ZeroAlign_By_Sex_", i,".csv"))
  Age0.perm[[i]] <- tmp
}

Age0.perm.df <- do.call(rbind, Age0.perm)
Age0.perm.df$diff <- Age0.perm.df$Median.M - Age0.perm.df$Median.F
Actual.age0$diff <- Actual.age0$Median.M - Actual.age0$Median.F

Pperm = sum(Age0.perm.df$diff > Actual.age0$diff) / 1000 # P=0.005

ggplot() +
  geom_histogram(aes(Age0.perm.df$diff, y = ..count..),
                 color = "black", fill = "#B4D3E7") +
  geom_vline(aes(xintercept = Actual.age0$diff), colour = "red", linetype="solid", linewidth=1)+
  labs(x = "Null distribution", y = "Frequency") +theme_classic()+
  #scale_y_discrete(breaks=NULL)+scale_x_continuous(breaks=NULL)+
  theme(axis.text=element_text(size=22, color="black"), 
        axis.title =element_text(size=22, color="black"),aspect.ratio = 1,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=22, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder,'/CV', CVthr,  '/Alignment_development/SexDifference_Permutation.tiff'), 
       width=15, height =14, units="cm")


