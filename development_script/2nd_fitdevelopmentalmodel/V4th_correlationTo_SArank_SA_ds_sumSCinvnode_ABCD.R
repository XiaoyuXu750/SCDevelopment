## Validation: large-scale matrix of 17 or 7
#### This script is to conduct correlation analysis 
#### between gam statistical indexes to connectional axis rank.
#### And draw scatter plots & matrix graphs. Supplementary Figure 5
library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)

rm(list = ls())
ds.resolution <- 7
elementnum <- ds.resolution*(ds.resolution+1) /2

demopath<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/demopath'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA', ds.resolution)

#### load data
CVthr = 75
gamresult<-readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_siteall_CV', CVthr,'_scale_TRUE.rds'))
gamresult$pfdr <- p.adjust(gamresult$bootstrap_pvalue, method="fdr")
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
perm.id.full<-readRDS(paste0("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA", ds.resolution, "_sphericalrotations_N10000.rds"))
#### source function
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))

#### 1. calculate correlation to SC rank
boxplot(gamresult$partialRsq)
gamresult <- within(gamresult, 
                    {partialRsq2 <- partialRsq
                    partialRsq2[which(partialRsq2>mean(partialRsq)+3*sd(partialRsq) | partialRsq2<mean(partialRsq)-3*sd(partialRsq))] <- NA
                    })
summary(gamresult$partialRsq2)
computevar <- "partialRsq2"
ds.resolution<-ds.resolution
SCrank_correlation <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE)
computevar <- "meanderv2"
SCrank_correlation <- rbind(SCrank_correlation, SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE))
# matrix of 17
#  ds.resolution Interest.var r.spearman  p.spin
# 1            17  partialRsq2 -0.6070170 0.00185
# 2            17    meanderv2  0.7507753 0.00130

# matrix  of 7
#   ds.resolution Interest.var r.spearman   p.spin
# 1             7  partialRsq2 -0.8030656 0.045225
# 2             7    meanderv2  0.8835363 0.003425

### 2. scatter plots
############################################
ds.resolution <- ds.resolution

## partial Rsq
computevar <- "partialRsq2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
summary(correlation.df$partialRsq2)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$partialRsq2
lwth = min(correlation.df$partialRsq2, na.rm = T); upth = max(correlation.df$partialRsq2, na.rm = T)
colorbarvalues.Rsq <- colorbarvalues(correlation.df$partialRsq2, (-lwth/(upth-lwth)))
rho <- round(SCrank_correlation$r.spearman[1], 2)
pspin <- round(SCrank_correlation$p.spin[1], 3)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq2, color=partialRsq2), size=5)+
  geom_smooth(aes(x=SCrank, y=partialRsq2), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.Rsq)+
  labs(x="Connectional axis rank", y=expression("Age effect (partial "*italic("R")^"2"*")"))+
  theme_classic()+
  theme(axis.text=element_text(size=24, color="black"), 
        axis.title =element_text(size=24),aspect.ratio = 1.05,axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.title = element_text(size=15, hjust = 0.5, vjust=0),
        plot.subtitle = element_text(size=21, hjust = 0.9, vjust=-6),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=17, height =14, units = "cm")

## mean 2nd derivatives
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
summary(gamresult$meanderv2)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2
lwth2 = min(correlation.df$meanderv2, na.rm = T); upth2 = max(correlation.df$meanderv2, na.rm = T)
colorbarvalues.meanderv2 <- colorbarvalues(correlation.df$meanderv2, (-lwth2/(upth2-lwth2)))
rho <- round(SCrank_correlation$r.spearman[2], 2)
pspin <- round(SCrank_correlation$p.spin[2], 3)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2, color=meanderv2), size=5)+
  geom_smooth(aes(x=SCrank, y=meanderv2), method ="lm", color="black", linewidth=1.2)+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.meanderv2)+
  labs(x="Connectional axis rank", y="Second derivative")+
  scale_y_continuous(breaks = c(-0.005,0,0.005,0.010), labels=c(-5,0,5,10))+
  theme_classic()+
  theme(axis.text=element_text(size=23.2, color="black"), 
        axis.title =element_text(size=23.2),aspect.ratio = 1.0,axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.title = element_text(size=15, hjust = 0.5, vjust=0),
        plot.subtitle = element_text(size=15, hjust = 0.1, vjust=-6),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_siteall.svg'), dpi=600, width=13, height =13, units = "cm")


### 3. matrix graphs for resolution at ds.resolution
############################################
Matrix.tmp <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)

linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))

computevarlist <- c("partialRsq2", "meanderv2")
colorprob <- c((-lwth/(upth-lwth)), (-lwth2/(upth2-lwth2)))
n=0
for (computevar in computevarlist){
  SCrank_correlation.df<-SCrankcorr(gamresult, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
  
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank_correlation.df[,2]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, ds.resolution)
  rownames(Matrix.tmp) <-seq(1, ds.resolution)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  n=n+1
  colorbarvalues.tmp <- colorbarvalues(matrixtmp.df.melt$value, colorprob[n])
  
  Matrix.tmp.sig <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
  Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (gamresult$pfdr<0.05)
  Matrix.tmp.sig[upper.tri(Matrix.tmp.sig)] <- t(Matrix.tmp.sig)[upper.tri(Matrix.tmp.sig)]
  colnames(Matrix.tmp.sig) <-seq(1, ds.resolution)
  rownames(Matrix.tmp.sig) <-seq(1, ds.resolution)
  matrixtmp.df.sig <- as.data.frame(Matrix.tmp.sig)
  matrixtmp.df.sig$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.sig.melt <- melt(matrixtmp.df.sig,id.vars=c("nodeid"))
  matrixtmp.df.sig.melt$variable<-as.numeric(matrixtmp.df.sig.melt$variable)
  matrixtmp.df.sig.melt$nodeid<-0-matrixtmp.df.sig.melt$nodeid
  matrixtmp.df.sig.melt$value<-as.numeric(matrixtmp.df.sig.melt$value)
  matrixtmp.df.sig.melt <- matrixtmp.df.sig.melt[-which(matrixtmp.df.sig.melt$value==0),]
  if (computevar=="partialRsq2"){
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
      scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
      geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=6)+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = ds.resolution+0.5 ,yend = -ds.resolution-0.5), color="black", linewidth=0.5)+
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
      scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
      scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "#053061")+
      #geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=6)+
      geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
      geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
      geom_segment(aes(x = 0.5 , y = -0.5 , xend = ds.resolution+0.5 ,yend = -ds.resolution-0.5), color="black", linewidth=0.5)+
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
  filename<-paste0(FigureFolder,'/CV', CVthr, "/Matrix", ds.resolution, "_sumSCinvnode_gamstats_Age8_22/", computevar, "_delLM_CV", CVthr,"_siteall.tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}









