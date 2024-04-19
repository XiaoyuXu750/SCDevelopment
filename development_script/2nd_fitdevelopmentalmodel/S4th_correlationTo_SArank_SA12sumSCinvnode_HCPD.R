#### This script is to conduct correlation analysis 
#### between gam statistical indexes to connectional axis rank.
#### And draw scatter plots & matrix graphs.
#### Figure 2. (a, c); Figure 3. (b, c); Supplementary Figure 2 (a,b); 
#### Supplementary Figure 4; Supplementary Figure 6 (a,b,d,e,g,h)
#### Spearman correlations were conducted.
#### The functions of spin tests came from Váša, F. et al. (2018, https://github.com/frantisekvasa/rotate_parcellation).
#### Network spin tests refer to Hansen et al., 2022, Nature Communications (https://github.com/netneurolab/hansen_crossdisorder_vulnerability/blob/main/code/03_disorder_similarity.py)
#### line 185 ~ 201
library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
rm(list = ls())
demopath<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/demopath'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_HCPD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'

#### load data
CVthr = 25
gamresult<-readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method="fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.merge.rds'))
perm.id.full<-readRDS("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA12_sphericalrotations_N10000.rds")
#### source function
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
#### description
meanSC <- colMeans(SCdata[,which(str_detect(names(SCdata), "SC."))])
meanlength <- colMeans(SCdata[,which(str_detect(names(SCdata), "length"))])
corr.test(meanSC, gamresult$partialRsq) # r=0.55, p=0
corr.test(meanSC, meanlength) # r=-0.78, p=0
#### convert critical ages of insignificantly developmental edges to NA
#### convert critical ages equal to age boundaries to NA
gamresult$meanlength <- meanlength
gamresult$increase.onset[gamresult$sig==FALSE]<-NA ; gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2,2)==8.08] <- NA
gamresult$increase.offset[gamresult$sig==FALSE]<-NA ; gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2,2)==21.92] <- NA
gamresult$peak.change[gamresult$sig==FALSE]<-NA
gamresult$peak.increase.change[gamresult$sig==FALSE]<-NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")
summary(gamresult)


#### 1. compute correlations to SC rank
computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change",
                     "meanderv2")
for (x in 1:5){
  computevar <- computevar.list[x]
  ds.resolution<-12
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE)
  if (x==1){
    SCrank_correlation <- correlation.df
  }else{
    SCrank_correlation <- rbind(SCrank_correlation, correlation.df)
  }
}
SCrank_correlation
# CV 75th
#     ds.resolution   Interest.var r.spearman  p.spin
# 1            12           partialRsq -0.1714174 0.21175
# 2            12      increase.onset2  0.2337374 0.073850
# 3            12     increase.offset2  0.4896016 0.012500
# 4            12 peak.increase.change  0.8661125 0.000375
# 5            12            meanderv2  0.8045298 0.002350

# CV 25th
#    ds.resolution   Interest.var r.spearman  p.spin
# 1            12           partialRsq -0.2806793 0.108300
# 2            12      increase.onset2  0.4309435 0.004650
# 3            12     increase.offset2  0.4932297 0.022125
# 4            12 peak.increase.change  0.8541169 0.000375
# 5            12            meanderv2  0.7882164 0.001375

#### 1.1 Validation: Control for fiber length
corr.test(gamresult$partialRsq, gamresult$meanlength, method = "pearson")
# r=-0.38,p=0
corr.test(gamresult$meanderv2, gamresult$meanlength, method = "pearson")
# r=0.53, p=0
meanderv2_control_length <- residuals(lm(meanderv2~meanlength, data=gamresult))
corr.test(meanderv2_control_length, gamresult$meanlength, method = "pearson") # r=0
gamresult$meanderv2_control_length <- meanderv2_control_length
SCrankcorr(gamresult, "meanderv2_control_length", 12, perm.id.full, dsdata=FALSE)
# CV75th
#   ds.resolution     Interest.var r.spearman  p.spin
# 1            12    meanderv2_control_length   0.768185 0.003725
# CV25th
#   ds.resolution     Interest.var r.spearman  p.spin
# 1            12    meanderv2_control_length   0.7721938 0.002125

partialRsq_control_length <- residuals(lm(partialRsq~meanlength, data=gamresult))
corr.test(partialRsq_control_length, gamresult$meanlength, method = "pearson") # r=0
gamresult$partialRsq_control_length <- partialRsq_control_length
SCrankcorr(gamresult, "partialRsq_control_length", 12, perm.id.full, dsdata=FALSE)
# CV75
#    ds.resolution      Interest.var r.spearman  p.spin
# 1            12    partialRsq_control_length -0.1444054 0.2077
# CV25
#    ds.resolution      Interest.var r.spearman  p.spin
# 1            12    partialRsq_control_length -0.2881279 0.142475

### 2. scatter plots
############################################
ds.resolution <- 12

## partial Rsq
computevar <- "partialRsq"

correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
#names(correlation.df) <- c("mat_a", "mat_b")
#write.csv(correlation.df, paste0(interfileFolder, '/correlationdf_onset.csv'), row.names = F)
summary(gamresult$partialRsq)
colorbarvalues.Rsq <- colorbarvalues(correlation.df$partialRsq, 0.5)
lmthr <- max(abs(gamresult$partialRsq))
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$partialRsq
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq, color=partialRsq, size=0.5))+
  geom_smooth(aes(x=SCrank, y=partialRsq), method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu",direction = -1, limit=c(-lmthr,lmthr))+  labs(x="Connectional axis rank", y="Partial R2")+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.8,
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder, '/CV', CVthr,'/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/CV', CVthr,'/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=17, height =14, units = "cm")

# increase onset
computevar <- "increase.onset2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$increase.onset2
colorbarvalues.onset <- colorbarvalues(correlation.df$increase.onset2, 0.5)
RdBu3 <- rev(brewer.pal(11, "RdBu"))
RdBu3 <- RdBu3[-c(1:3)]
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=increase.onset2, color=increase.onset2, size=0.5))+
  geom_smooth(aes(x=SCrank, y=increase.onset2), method ="lm", color="black")+
  #scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.onset)+
  scale_colour_gradientn(colours = RdBu3,values=colorbarvalues.onset, space="Lab")+
  labs(x="Connectional axis rank", y="Onset age (years)")+
  scale_x_continuous(limits = c(15,78))+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.8,
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=17, height =14, units = "cm")

# increase offset
computevar <- "increase.offset2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$increase.offset2
colorbarvalues.offset <- colorbarvalues(correlation.df$increase.offset2, 0.5)
RdBu4 <- rev(brewer.pal(11, "RdBu"))
RdBu4 <- RdBu4[-c(9:11)]
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=increase.offset2, color=increase.offset2, size=0.5))+
  geom_smooth(aes(x=SCrank, y=increase.offset2), method ="lm", color="black")+
  #scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.offset)+
  scale_colour_gradientn(colours = RdBu4,values=colorbarvalues.offset, space="Lab")+
  labs(x="Connectional axis rank", y="Offset age (years)")+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.8,
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=17, height =14, units = "cm")

# peak age
computevar <- "peak.increase.change"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$peak.increase.change
colorbarvalues.onset <- colorbarvalues(correlation.df$peak.increase.change, 0.5)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=peak.increase.change, color=peak.increase.change, size=0.5))+
  geom_smooth(aes(x=SCrank, y=peak.increase.change), method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.onset)+
  labs(x="Connectional axis rank", y="Age of max slope (years)")+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20),aspect.ratio = 0.8,
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=17, height =14, units = "cm")

# meanderv2
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2
colorbarvalues.meanderv2 <- colorbarvalues(correlation.df$meanderv2, 0.5)
if (CVthr==75){
  mytheme <- theme(axis.text=element_text(size=20, color="black"), 
                   axis.title =element_text(size=20),aspect.ratio = 0.75,
                   axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
                   plot.title = element_text(size=20, hjust = 0.5, vjust=2),
                   plot.background=element_rect(fill="transparent"),
                   panel.background=element_rect(fill="transparent"),
                   legend.position = "none")
  width =17.5; height = 15
}else{
  mytheme <- theme(axis.text=element_text(size=23, color="black"), 
                   axis.title =element_text(size=23),aspect.ratio = 1,
                   axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
                   plot.title = element_text(size=20, hjust = 0.5, vjust=2),
                   plot.background=element_rect(fill="transparent"),
                   panel.background=element_rect(fill="transparent"),
                   legend.position = "none")
  width = height = 14
}

ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2, color=SCrank), size=5)+
  geom_smooth(aes(x=SCrank, y=meanderv2),linewidth=2, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="Connectional axis rank", y="Second derivative")+
  scale_y_continuous(breaks = c(-0.003, 0, 0.003), labels = c(-3, 0, 3))+
  theme_classic()+mytheme
  
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=width, height =height, units = "cm")

# meanderv2_controllength
computevar <- "meanderv2_control_length"
correlation.df <- SCrankcorr(gamresult,computevar, ds.resolution, perm.id.full,dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2_control_length
colorbarvalues.meanderv2_control_length <- colorbarvalues(correlation.df$meanderv2_control_length, 0.5)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2_control_length, color=meanderv2_control_length, size=0.2))+
  geom_smooth(aes(x=SCrank, y=meanderv2_control_length),linewidth=1.2, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu",values =colorbarvalues.meanderv2_control_length, direction = -1)+
  labs(x="Connectional axis rank", y="Second derivative")+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=21, color="black"), 
        axis.title =element_text(size=21),aspect.ratio = 0.9,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=15, height =14, units = "cm")


### 3. matrix graphs for resolution of 12
############################################
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)
computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2", "meanderv2_control_length")
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))
SCrank_correlation.df<-mclapply(1:6, function(x){
  computevar <- computevar.list[x]
  ds.resolution<-12
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
  return(correlation.df)
}, mc.cores = 4)
colorbar.prob <- c(0.5, 0.4, 0.6, 0.5, 0.5, 0.5)

for (i in 1:6){
  computevar <- computevar.list[i]
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank_correlation.df[[i]][,2]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, 12)
  rownames(Matrix.tmp) <-seq(1, 12)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, 12)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  colorbarvalues.tmp <- colorbarvalues(SCrank_correlation.df[[i]][,2], colorbar.prob[i])
  if (computevar=="partialRsq"){
    lmthr <- max(abs(gamresult$partialRsq))
    Matrix.tmp.sig <- matrix(NA, nrow = 12, ncol=12)
    Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (gamresult$pfdr<0.05)
    Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
    colnames(Matrix.tmp.sig) <-seq(1, 12)
    rownames(Matrix.tmp.sig) <-seq(1, 12)
    matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
    matrixtmp.df.sig$nodeid <- seq(1, 12)
    matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
    matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
    matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
    matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
    matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value==0),]
    
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
      scale_color_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
      geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=6)+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
      ggtitle(label = computevar)+labs(x=NULL, y=NULL)+
      scale_y_continuous(breaks=NULL, labels = NULL)+
      scale_x_continuous(breaks=NULL, labels = NULL)+
      theme(axis.line = element_blank(), 
            #axis.ticks=element_line(linewidth = 0),
            axis.text.x=element_text(size=12, angle=45, hjust=1), 
            axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
            axis.title =element_text(size=18),
            plot.title = element_text(size=18, hjust = 0.5),
            legend.title=element_text(size=18),
            legend.text=element_text(size=18), 
            panel.background=element_rect(fill=NA),
            panel.grid.major=element_line(linewidth = 0), 
            panel.grid.minor=element_line(linewidth = 1))
    
  }else{
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      #paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.tmp,oob = squish) +
      #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1,values=colorbarvalues.tmp, oob = squish) +
      scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "grey")+
      scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "grey")+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
      ggtitle(label = computevar)+labs(x=NULL, y=NULL)+
      scale_y_continuous(breaks=NULL, labels = NULL)+
      scale_x_continuous(breaks=NULL, labels = NULL)+
      theme(axis.line = element_blank(), 
            #axis.ticks=element_line(linewidth = 0),
            axis.text.x=element_text(size=12, angle=45, hjust=1), 
            axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
            axis.title =element_text(size=18),
            plot.title = element_text(size=18, hjust = 0.5),
            legend.title=element_text(size=18),
            legend.text=element_text(size=18), 
            panel.background=element_rect(fill=NA),
            panel.grid.major=element_line(linewidth = 0), 
            panel.grid.minor=element_line(linewidth = 1))
  }
  
  Fig
  filename<-paste0(FigureFolder,'/CV', CVthr,  "/Matrix12_sumSCinvnode_gamstats_Age8_22/", computevar, "_12net_delLM_CV75.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}



