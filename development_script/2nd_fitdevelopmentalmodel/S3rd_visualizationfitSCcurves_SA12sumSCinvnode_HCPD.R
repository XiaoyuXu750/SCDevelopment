## This script is to generate fitted values from gam models.
## The predicted SC strength will be generated as the age was sampled from 8 to 22, 
## with a total of 1000 data points, covariables will be set as median or mode.
## The data will be used to draw developmental trajectories. (Figure 2(b, d, e, f); 
## Figure 3(d); Supplementary Figure 2 (c))
library(R.matlab)
library(mgcv)
library(psych)
library(tidyverse)
library(parallel)
library(scales)
library(openxlsx)
library(gratia)
library(RColorBrewer)
library(paletteer)
rm(list = ls())

resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_HCPD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'

#### load data
CVthr = 75
gamresultsum.SAorder.delLM<-readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/plotdata_derivatives.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
#### generate fitted values for developmental trajectories
plotdatasum<-mclapply(1:78, function(x){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, "age")
  return(plotdata)
}, mc.cores = 2)
plotdatasum.df<-as.data.frame(matrix(NA, nrow = 1, ncol=19))
names(plotdatasum.df)<-c(names(plotdatasum[[2]])[1:15],"SC_label",  "SCrank", 
                         "PartialRsq", "meanderiv2")
#### SA12 index & SC rank
Matrix12<-matrix(NA, nrow=12, ncol=12)
indexup12 <- upper.tri(Matrix12)
indexsave12 <- !indexup12
Matrix12.index<-Matrix12
#13*12/2=78
Matrix12.index[indexsave12]<-c(1:78)
#S-A connectional axis rank
Matrix12.SCrank<-Matrix12
for (x in 1:12){
  for (y in 1:12){
    Matrix12.SCrank[x,y]<-x^2+y^2
  }
}
Matrix12.SCrank[indexup12]<-NA
Matrix12.SCrank[indexsave12]<-rank(Matrix12.SCrank[indexsave12], ties.method = "average")
# rbind plotdata
for (i in 1:78){
  tmp<-plotdatasum[[i]][,-16]
  tmp$SC_label<-names(plotdatasum[[i]])[16]
  tmp$SCrank<-Matrix12.SCrank[indexsave12][i]
  tmp$PartialRsq<-gamresultsum.SAorder.delLM$partialRsq[i]
  tmp$meanderiv2<-gamresultsum.SAorder.delLM$meanderv2[i]
  plotdatasum.df<-rbind(plotdatasum.df, tmp)
}
plotdatasum.df<-plotdatasum.df[-1, ]

## plot
# Figure 2 (b)
lmthr <- max(abs(gamresultsum.SAorder.delLM$partialRsq))
ggplot()+
  geom_line(data=plotdatasum.df, aes(x=age, y=fit.ratio, group=SC_label, color=PartialRsq), size=0.8, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1, limit=c(-lmthr,lmthr))+
  labs(x="Age (years)", y="SC strength (ratio)")+
  theme_classic()+
  theme(axis.text=element_text(size=20.5, color="black"), 
        axis.title =element_text(size=20.5, color="black"),aspect.ratio = 1,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_sumSCinvnode_fit/devcurve_Rsq_fit.ratio.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_sumSCinvnode_fit/devcurve_Rsq_fit.ratio.svg'), dpi=600, width=15, height =15, units = "cm")
# Figure 2 (d)
colorbarvalues.meanderiv2 <- colorbarvalues(plotdatasum.df$meanderiv2, 0.5)
SC_label_derv2_order <- gamresultsum.SAorder.delLM$parcel[order(gamresultsum.SAorder.delLM$meanderv2)]
plotdatasum.df$SC_label2 <- factor(plotdatasum.df$SC_label, levels=SC_label_derv2_order)
ggplot()+
  geom_line(data=plotdatasum.df, aes(x=age, y=fit.Z, group=SC_label2, color=meanderiv2), size=0.8, alpha=0.8)+
  #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.meanderiv2, oob = squish) +
  scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.meanderiv2, direction = -1)+
  labs(x="Age (years)", y="SC strength (z-score)")+
  #scale_color_manual(values = rev(brewer.pal(10, "RdBu")))+
  scale_y_continuous(breaks = c(-1.5, 0.0, 1.5))+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23, color="black"),aspect.ratio = 0.97,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_sumSCinvnode_fit/devcurve_meanderv2_fit.Z.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_sumSCinvnode_fit/devcurve_meanderv2_fit.Z.svg'), dpi=600, width=16, height =15, units = "cm")

#### Example, Figure 2 (e, f)
############
BuRd <- rev(brewer.pal(10, "RdBu"))
##visualization
# i=75
# i=36
# i=2
i=77
SClabel<-paste0("SC.",i, "_h")
range <- max(gamresultsum.SAorder.delLM$meanderv2) - min(gamresultsum.SAorder.delLM$meanderv2)
coloridx=round((gamresultsum.SAorder.delLM$meanderv2[gamresultsum.SAorder.delLM$parcel==SClabel]-min(gamresultsum.SAorder.delLM$meanderv2))/range * 10)
if (coloridx==0){coloridx=1}
colorID <- BuRd[coloridx]

plotdatasum.df.tmp<-plotdatasum.df[plotdatasum.df$SC_label==SClabel,]
scatterdata<-gammodelsum[[i]]$model
scatterdata$SC <- scatterdata[,1]
Scatter_Fig<-ggplot()+
  geom_ribbon(data=plotdatasum.df.tmp, aes(x=age, ymin=selo, ymax=sehi), alpha=0.3, fill=colorID)+
  geom_point(data=scatterdata, aes(x=age, y=SC), alpha=0.5, color=colorID)+
  geom_line(data=plotdatasum.df.tmp, aes(x=age, y=fit), size=1.4, alpha=1, color=colorID)+
  labs(y="SC strength (ratio)")+xlab(NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=19, color="black"), 
        axis.title =element_text(size=19),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        aspect.ratio = 0.8,axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20))
Scatter_Fig
# derivative plot
deriv.SA12.tmp<-derivative[derivative$label_ID==SClabel, ]
deriv.SA12.tmp$h <- 1
derivplot <-ggplot(data=deriv.SA12.tmp)+
  geom_bar(aes(x=age, y=1, fill = significant.derivative, color=significant.derivative),stat = "identity", position = "stack")+
  scale_fill_gradient2(high = colorID, low = "white",midpoint=0, na.value = "white", labels=NULL) +
  scale_color_gradient2(high = colorID, low = "white",midpoint=0, na.value = "white",  labels=NULL) +
  scale_y_continuous(breaks = NULL)+
  ylab(NULL)+xlab(NULL)+
  scale_x_continuous(breaks = NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=20, color='black'),
        axis.title = element_text(size = 20),
        axis.line.y=element_line(size=0),
        axis.line.x=element_line(size=0),
        legend.position = 'none')

derivplot
# merge plots
allplots <- list(Scatter_Fig,derivplot)
mergeplot <- cowplot::plot_grid(rel_heights = c(16,1), plotlist = allplots, 
                                align = "v", axis = "lr",greedy=F, ncol = 1, nrow = 2,
                                hjust = -0.5)
mergeplot

ggsave(paste0(FigureFolder,'/CV',CVthr, '/SA12_sumSCinvnode_fit/SA12_delLM_', SClabel, '_CV75.tiff'),mergeplot, width=14, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV',CVthr, '/SA12_sumSCinvnode_fit/SA12_delLM_', SClabel, '_CV75.svg'),mergeplot, dpi = 600, width=14, height =14
       , units = "cm")
#####################

## Average fitted values for 10 deciles of connectional axis
# Figure 3 (d)
SA12_10 <- data.frame(SCrank=Matrix12.SCrank[indexsave12])
SA12_10 <- SA12_10 %>%
  mutate(decile = ntile(SCrank, 10))
table(SA12_10$decile)
SA12_10$SC_label <- paste0("SC.", c(1:78), "_h")
write.csv(SA12_10, paste0(interfileFolder, '/SA12_10.csv'), row.names = F)
plotdatasum.df.label <- merge(plotdatasum.df, SA12_10, by="SC_label", all.x=T)
names(plotdatasum.df.label)
plotdatasum.df.decile <- plotdatasum.df.label %>%
  group_by(decile, age) %>%
  summarise(fit.avg = mean(fit), SCranktype_order=mean(decile))

names(plotdatasum.df.decile)
plotdatasum.df.decile <- plotdatasum.df.decile %>%
  group_by(decile) %>%
  mutate(fit.Z = scale(fit.avg), fit.ratio=fit.avg/fit.avg[1])

ggplot(data=plotdatasum.df.decile, aes(x=age, y=fit.Z, group=decile, color=decile))+
  geom_line(size=1.5, alpha=0.8)+
  #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.meanderiv2, oob = squish) +
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="Age (years)", y="SC strength (z-score)")+
  #scale_color_manual(values = rev(brewer.pal(6, "RdBu")))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"),aspect.ratio = 0.8,
        plot.background=element_rect(fill="transparent"),legend.position = "none",
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5))
ggsave(paste0(FigureFolder,'/CV',CVthr,  '/SA12_decile_sumSCinvnode_fit/devcurve_SCrank_fit.Z_SCtype10.tiff'), dpi=600, width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV',CVthr,  '/SA12_decile_sumSCinvnode_fit/devcurve_SCrank_fit.Z_SCtype10.svg'), dpi=600, width=15, height =13, units = "cm")




