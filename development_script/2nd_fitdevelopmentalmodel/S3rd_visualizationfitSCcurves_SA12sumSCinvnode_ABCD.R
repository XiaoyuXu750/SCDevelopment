## This script is to generate fitted values from gam models.
## The predicted SC strength will be generated as the age was sampled from 8.9 to 13.8, 
## with a total of 1000 data points, covariables will be set as median or mode.
## The data will be used to draw developmental trajectories. (Figure 5(b, e), Sup-Figure 3(c))
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

resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA12'

#### load data
CVthr = 75
gamresultsum.SAorder.delLM<-readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
gammodelsum<-readRDS(paste0(interfileFolder, '/gammodel78_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
#### source function
source(paste0(functionFolder, '/plotdata_generate.R'))
source(paste0(functionFolder, '/plotdata_derivatives.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
source(paste0(functionFolder, '/gammsmooth.R'))
#### generate function data
#### SA12 index & SC rank
Matrix12<-matrix(NA, nrow=12, ncol=12)
indexup12 <- upper.tri(Matrix12)
indexsave12 <- !indexup12
Matrix12.index<-Matrix12
#13*12/2=78
Matrix12.index[indexsave12]<-c(1:78)
#SC rank
Matrix12.SCrank<-Matrix12
for (x in 1:12){
  for (y in 1:12){
    Matrix12.SCrank[x,y]<-(x+y)^2+(x-y)^2
  }
}
Matrix12.SCrank[indexup12]<-NA
Matrix12.SCrank[indexsave12]<-rank(Matrix12.SCrank[indexsave12], ties.method = "average")
gamresultsum.SAorder.delLM$SCrank <- Matrix12.SCrank[indexsave12]
# plot data
plotdatasum<-mclapply(1:78, function(x){
  modobj<-gammodelsum[[x]]
  plotdata<- plotdata_generate(modobj, "age")
  plotdata$SC_label <- names(plotdata)[16]
  plotdata[,16] <- NULL
  plotdata$SCrank <- Matrix12.SCrank[indexsave12][x]
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
# Figure 5 (b)
colorbarvalues.Rsq <- colorbarvalues(plotdatasum.df$partialRsq2, 0.4)
RColorBrewer::brewer.pal(11, "RdBu")
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

ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_sumSCinvnode_fit/devcurve_SCrank_fit.ratio.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_sumSCinvnode_fit/devcurve_Rsq_fit.ratio.svg'), dpi=600, width=14, height =13, units = "cm")

#### Example
############
BuRd <- rev(brewer.pal(10, "RdBu"))
##visualization
# i=77
# i=3
i=3
SClabel<-paste0("SC.",i, "_h")
gamresultsum.SAorder.delLM$partialRsq2 <- gamresultsum.SAorder.delLM$partialRsq
gamresultsum.SAorder.delLM$partialRsq2[which.min(gamresultsum.SAorder.delLM$partialRsq2)] <- NA
range <- range.default(gamresultsum.SAorder.delLM$partialRsq2, na.rm = T)[2] - range.default(gamresultsum.SAorder.delLM$partialRsq2, na.rm = T)[1]
coloridx=round(gamresultsum.SAorder.delLM$SCrank[gamresultsum.SAorder.delLM$parcel==SClabel]/78 * 10)
coloridx=9
if (coloridx==0){coloridx=1}
colorID <- BuRd[coloridx]
plotdatasum.df.tmp<-plotdatasum.df[plotdatasum.df$SC_label==SClabel,]
plotdatasum.df.tmp$selo.ratio <- plotdatasum.df.tmp$selo / plotdatasum.df.tmp$fit[1]
plotdatasum.df.tmp$sehi.ratio <- plotdatasum.df.tmp$sehi / plotdatasum.df.tmp$fit[1]
scatterdata<-gammodelsum[[i]]$gam$model
scatterdata$SC <- scatterdata[,1] / plotdatasum.df.tmp$fit[1]
model.tmp <- gammodelsum[[i]]
model.tmp$gam$data <- scatterdata
pointdata <- visreg(model.tmp$gam, 'age', plot=F)
scatterdata$SCpoint <- pointdata$res$visregRes

Scatter_Fig<-ggplot()+
  geom_ribbon(data=plotdatasum.df.tmp, aes(x=age, ymin=selo, ymax=sehi), alpha=0.3, fill=colorID)+
  geom_line(data =scatterdata, aes(x = age,y =SCpoint, group = subID),alpha = .2,color="grey", size = 0.5)+
  geom_jitter(data=scatterdata, aes(x=age, y=SCpoint),color="black",shape=16, size=1, alpha=0.3)+
  geom_line(data=plotdatasum.df.tmp, aes(x=age, y=fit), size=1.4, alpha=1, color=colorID)+
  scale_y_continuous(limits = c(0.4, 2))+
  labs(y="SC strength (ratio)")+xlab(NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=24, color="black"), 
        axis.title =element_text(size=24),
        plot.title = element_text(size=24, hjust = 0.5, vjust=2),
        aspect.ratio = 1,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20))
Scatter_Fig
# derivative plot
deriv.SA12.tmp<-derivative[derivative$label_ID==SClabel, ]
deriv.SA12.tmp$h <- 1

derivplot<-ggplot(data=deriv.SA12.tmp)+
  geom_bar(aes(x=age, y=1, fill = significant.derivative, color=significant.derivative),stat = "identity", position = "stack")+
  scale_fill_gradient2(low ="white" , high  = colorID,midpoint=0, na.value = "white", labels=NULL) +
  scale_color_gradient2(low ="white" , high  = colorID,midpoint=0, na.value = "white",  labels=NULL) +
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
ggsave(paste0(FigureFolder,'/CV',CVthr, '/SA12_sumSCinvnode_fit/SA12_delLM_', SClabel, '_CV75.svg'),mergeplot, dpi = 600, width=15, height =15, units = "cm")
#####################

# decile
# Figure 5 (e)
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
plotdatasum.df.decile <- plotdatasum.df.decile %>%
  group_by(decile) %>%
  mutate(fit.Z = scale(fit.avg), fit.ratio=fit.avg/fit.avg[1])
ggplot(data=plotdatasum.df.decile, aes(x=age, y=fit.Z, group=decile, color=decile))+
  geom_line(size=1.5, alpha=0.8)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="Age (years)", y="SC strength (ratio)")+
  #scale_y_continuous(limits = c(0.94, 1.095), breaks=c(1.0,1.1))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"),aspect.ratio = 0.85,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=20, hjust = 0.5), legend.position = "none")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_decile_sumSCinvnode_fit/devcurve_SCrank_fit.Z_SCtype.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV',CVthr, '/SA12_decile_sumSCinvnode_fit/devcurve_SCrank_fit.Z_SCtype.svg'), dpi=600, width=14, height =13, units = "cm")

mat.tmp <- matrix(c(1:10), 1, 10)
image(mat.tmp, col=brewer.pal(11, "RdBu")[2:11])

