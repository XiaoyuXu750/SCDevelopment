## This script is to generate fitted values from gam models.
## The predicted SC strength will be generated as the age was sampled from 8.9 to 13.8, 
## with a total of 1000 data points, covariables will be set as median or mode.
## The data will be used to draw developmental trajectories. (Supplementary Figure 5)
library(R.matlab)
library(mgcv)
library(psych)
library(tidyverse)
library(parallel)
library(scales)
library(gratia)
library(RColorBrewer)
library(visreg)
rm(list = ls())
ds.resolution <- 17
elementnum <- ds.resolution*(ds.resolution+1) /2

resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA', ds.resolution)

#### load data
CVthr = 75
gamresultsum.SAorder.delLM<-readRDS(paste0(interfileFolder, '/gamresults',elementnum,'_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel',elementnum,'_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution,'_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df',elementnum,'_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/plotdata_derivatives.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
source(paste0(functionFolder, '/gammsmooth.R'))
#### generate function data
#### SAds.resolution index & SC rank
Matrixds.resolution<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
indexupds.resolution <- upper.tri(Matrixds.resolution)
indexsaveds.resolution <- !indexupds.resolution
Matrixds.resolution.index<-Matrixds.resolution
#(ds.resolution+1)*ds.resolution/2=elementnum
Matrixds.resolution.index[indexsaveds.resolution]<-c(1:elementnum)
#SC rank
Matrixds.resolution.SCrank<-Matrixds.resolution
for (x in 1:ds.resolution){
  for (y in 1:ds.resolution){
    Matrixds.resolution.SCrank[x,y]<-(x+y)^2+(x-y)^2
  }
}
Matrixds.resolution.SCrank[indexupds.resolution]<-NA
Matrixds.resolution.SCrank[indexsaveds.resolution]<-rank(Matrixds.resolution.SCrank[indexsaveds.resolution], ties.method = "average")
gamresultsum.SAorder.delLM$SCrank <- Matrixds.resolution.SCrank[indexsaveds.resolution]
# plot data
plotdatasum<-mclapply(1:elementnum, function(x){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, "age")
  plotdata$SC_label <- names(plotdata)[16]
  plotdata[,16] <- NULL
  plotdata$SCrank <- Matrixds.resolution.SCrank[indexsaveds.resolution][x]
  plotdata$PartialRsq <- gamresultsum.SAorder.delLM$partialRsq[x]
  plotdata$meanderv2 <- gamresultsum.SAorder.delLM$meanderv2[x]
  return(plotdata)
}, mc.cores = 2)
plotdatasum.df <- do.call(rbind, lapply(plotdatasum, as.data.frame))
summary(plotdatasum.df$age)
plotdatasum.df$partialRsq2 <- plotdatasum.df$PartialRsq
min(gamresultsum.SAorder.delLM$partialRsq)
plotdatasum.df$partialRsq2[which(plotdatasum.df$partialRsq2==min(plotdatasum.df$partialRsq2, na.rm = T))] <- NA
summary(plotdatasum.df$partialRsq2)
## plot
# Supplemetary Figure 5
lwth = min(plotdatasum.df$partialRsq2, na.rm = T); upth = max(plotdatasum.df$partialRsq2, na.rm = T)
colorbarvalues.Rsq <- colorbarvalues(plotdatasum.df$partialRsq2, (-lwth/(upth-lwth)))

ggplot()+
  geom_line(data=plotdatasum.df, aes(x=age, y=fit.ratio, group=SC_label, color=partialRsq2), size=0.8, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.Rsq,na.value="#053061", direction = -1)+
  labs(x="Age (years)", y="SC strength (ratio)")+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"),aspect.ratio = 0.85,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=20, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA', ds.resolution,'_sumSCinvnode_fit/devcurve_SCrank_fit.ratio.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA', ds.resolution,'_sumSCinvnode_fit/devcurve_Rsq_fit.ratio.svg'),  dpi=600, width=14, height =13, units = "cm")

# decile
# Supplemetary Figure 5
SAds.resolution_10 <- data.frame(SCrank=Matrixds.resolution.SCrank[indexsaveds.resolution])
SAds.resolution_10 <- SAds.resolution_10 %>%
  mutate(decile = ntile(SCrank, 10))
table(SAds.resolution_10$decile)
SAds.resolution_10$SC_label <- paste0("SC.", c(1:elementnum), "_h")
plotdatasum.df.label <- merge(plotdatasum.df, SAds.resolution_10, by="SC_label", all.x=T)
names(plotdatasum.df.label)
plotdatasum.df.decile <- plotdatasum.df.label %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(fit), SCranktype_order=mean(decile))
ggplot(data=plotdatasum.df.decile, aes(x=age, y=fit.avg, group=decile, color=decile))+
  geom_line(size=1.5, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="Age (years)", y="SC strength (ratio)")+
  scale_y_continuous(limits = c(0.9, 1.095), breaks=c(0.9, 1.0,1.1))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"),aspect.ratio = 0.9,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=20, hjust = 0.5), legend.position = "none")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA', ds.resolution,'_decile_sumSCinvnode_fit/devcurve_SCrank_fit.ratio_SCtype.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA', ds.resolution,'_decile_sumSCinvnode_fit/devcurve_SCrank_fit.ratio_SCtype.svg'), dpi=600, width=13, height =12, units = "cm")

mat.tmp <- matrix(c(1:10), 1, 10)
image(mat.tmp, col=brewer.pal(11, "RdBu")[2:11])
