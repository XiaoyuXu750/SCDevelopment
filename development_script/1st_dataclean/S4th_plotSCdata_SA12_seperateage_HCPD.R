## This script is to plot SC matrix at specific age
library(R.matlab);
library(ggplot2);
library(tidyverse)
library(corrplot)
library(reshape)
rm(list = ls())

demopath<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/demopath'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCP'
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCP_final/SA12'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'

schaefer400_index_SA<-read.csv(paste0(interfileFolder, '/schaefer400_index_SA.SAorder.delLM.csv'))

# load data
CVthr = 25
SCdata.sum.merge.over8 <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.merge.rds'))
source(paste0(functionFolder, '/colorbarvalue.R'))
age.time <- c(8, 15, 21)
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))
meanSCdata.sepage <- list()
for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum.merge.over8$age>=age.tmp & SCdata.sum.merge.over8$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum.merge.over8[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:79)])
  meanSCdata.sepage[[i]] <- meanSCdata.tmp}
lwth<-min(c(meanSCdata.sepage[[1]], meanSCdata.sepage[[2]], meanSCdata.sepage[[3]]))#0.411828
upth<-max(c(meanSCdata.sepage[[1]], meanSCdata.sepage[[2]], meanSCdata.sepage[[3]]))#12.60254
for (i in 1:3){
  age.tmp <- age.time[i]
  index.tmp <- which(SCdata.sum.merge.over8$age>=age.tmp & SCdata.sum.merge.over8$age<(age.tmp+1))
  SCdata.tmp <- SCdata.sum.merge.over8[index.tmp, ]
  meanSCdata.tmp <- colMeans(SCdata.tmp[,c(2:79)])
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- meanSCdata.tmp
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, 12)
  rownames(Matrix.tmp) <-seq(1, 12)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, 12)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  Fig<-ggplot(data =matrixtmp.df.melt)+
    geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
    scale_fill_distiller(type="seq", palette = "RdBu",limit=c(lwth, upth), na.value = "grey")+
    scale_color_distiller(type="seq", palette = "RdBu",limit=c(lwth, upth),na.value = "grey")+
    geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
    geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
    geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
    ggtitle(label = paste("age", age.tmp, "~", (age.tmp+1)))+labs(x=NULL, y=NULL)+
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
  Fig
  filename<-paste0(FigureFolder, "/Matrix12_sumSCinvnode_Age8_22/Age_", age.tmp, "_12net_delLM_CV", CVthr,".tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}

difmeanSC_15_8 <- (meanSCdata.sepage[[2]] - meanSCdata.sepage[[1]]) / meanSCdata.sepage[[1]]
difmeanSC_21_15 <- (meanSCdata.sepage[[3]] - meanSCdata.sepage[[2]]) / meanSCdata.sepage[[2]]
lwth <- min(c(difmeanSC_15_8, difmeanSC_21_15)); upth <- max(c(difmeanSC_15_8, difmeanSC_21_15))
## Mean SC 15 - Mean SC 8
Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- difmeanSC_15_8
Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
colnames(Matrix.tmp) <-seq(1, 12)
rownames(Matrix.tmp) <-seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
values1 <- colorbarvalues(matrixtmp.df.melt$value, 0.9)
Fig<-ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
  scale_fill_distiller(type="seq", palette = "RdBu",na.value = "grey", values=values1)+
  scale_color_distiller(type="seq", palette = "RdBu",na.value = "grey", values=values1)+
  geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
  geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
  geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
  ggtitle(label = "SC at age 15 minus age 8")+labs(x=NULL, y=NULL)+
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
Fig
filename<-paste0(FigureFolder, "/Matrix12_sumSCinvnode_Age8_22/Age_15_minus_8_12net_delLM_CV75_perct.tiff")
ggsave(filename, Fig,  height = 18, width = 20, units = "cm")

## Mean SC 21 - Mean SC 15
Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- difmeanSC_21_15
Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
colnames(Matrix.tmp) <-seq(1, 12)
rownames(Matrix.tmp) <-seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
values2 <- colorbarvalues(matrixtmp.df.melt$value, 0.5)
Fig<-ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
  scale_fill_distiller(type="seq", palette = "RdBu", values=values2,na.value = "grey")+
  scale_color_distiller(type="seq", palette = "RdBu",values=values2,na.value = "grey")+
  geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
  geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
  geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
  ggtitle(label = "SC at age 21 minus age 15")+labs(x=NULL, y=NULL)+
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
Fig
filename<-paste0(FigureFolder, "/Matrix12_sumSCinvnode_Age8_22/Age_22_minus_15_12net_delLM_CV75_perct.tiff")
ggsave(filename, Fig,  height = 18, width = 20, units = "cm")

## mean SC with 12 parcel bin indicating nodes
meanSCdata.tmp <- colMeans(SCdata.sum.merge.over8[,c(2:79)])
Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- meanSCdata.tmp
Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
colnames(Matrix.tmp) <-seq(1, 12)
rownames(Matrix.tmp) <-seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)

Fig<-ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
  scale_fill_distiller(type="seq", palette = "RdBu",limit=c(0.4, 13), na.value = "grey")+
  scale_color_distiller(type="seq", palette = "RdBu",limit=c(0.4, 13),na.value = "grey")+
  geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
  geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
  geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
  labs(x=NULL, y=NULL)+
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
Fig
filename<-paste0(FigureFolder, "/Matrix12_sumSCinvnode_Age8_22/Age_all_12net_delLM_CV75.tiff")
ggsave(filename, Fig,  height = 18, width = 20, units = "cm")

# 12*12 --> 78
meanSCdata.tmp <- meanSCdata.sepage[[3]]
ggplot()+
  geom_tile(aes(x=1, y=c(-1:-78), fill = meanSCdata.tmp, color=meanSCdata.tmp))+
  scale_fill_distiller(type="seq", palette = "RdBu",limit=c(0.4, 13),na.value = "grey")+
  scale_color_distiller(type="seq", palette = "RdBu",limit=c(0.4, 13),na.value = "grey")+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(breaks=NULL, labels = NULL)+
  scale_x_continuous(breaks=NULL, labels = NULL)+
  theme(axis.line = element_blank(), 
        #axis.ticks=element_line(linewidth = 0),
        axis.text.x=element_text(size=12, angle=45, hjust=1), 
        axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
        axis.title =element_text(size=18),
        plot.title = element_text(size=18, hjust = 0.5),
        legend.position="none",
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(linewidth = 0), 
        panel.grid.minor=element_line(linewidth = 1))
ggsave(paste0(FigureFolder, "/SC78_bin_age21.tiff"), width=1, height = 12, units = "cm")


# bin
schaefer400_index_SA$network_label <- factor(schaefer400_index_SA$network_label, levels=c("Vis",
                                                                                  "SomMot", "DorsAttn", "SalVentAttn","Limbic", "Cont", "Default"))
schaefer400_index_SA$nodeid = c(-1:-376)
tab12 <-table(schaefer400_index_SA$SA12node)
sepid <- cumsum(tab12)
linerange_frame1<-data.frame(y=c(-0.5, -sepid[1:11]+0.5, -sepid[12]-0.5), xmin =rep(0.5, times=13), xmax =rep(1.5, times=13))
linerange_frame2<-data.frame(x=c(0.5,1.5), ymin=c(-376.5, -376.5), ymax=c(-0.5, -0.5))

ggplot()+
  geom_tile(data =schaefer400_index_SA, aes(x=1, y=nodeid, fill = network_label, color=network_label))+
  scale_fill_manual(values=c("#823D8A", "#6E9EC3", "#2E7A3D", "#D867FF", "grey50", "#DAA34B",
                             "#D76779"))+
  scale_color_manual(values=c("#823D8A", "#6E9EC3", "#2E7A3D", "#D867FF", "grey50", "#DAA34B",
                              "#D76779"))+
  geom_linerange(data=linerange_frame1, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
  geom_linerange(data=linerange_frame2, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(breaks=NULL, labels = NULL)+
  scale_x_continuous(breaks=NULL, labels = NULL)+
  theme(axis.line = element_blank(), 
        #axis.ticks=element_line(linewidth = 0),
        axis.text.x=element_text(size=12, angle=45, hjust=1), 
        axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
        axis.title =element_text(size=18),
        plot.title = element_text(size=18, hjust = 0.5),
        legend.position="none",
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(linewidth = 0), 
        panel.grid.minor=element_line(linewidth = 1))
ggsave(paste0(FigureFolder, "/SA12_Yeo7_bin.tiff"), width=1, height = 12, units = "cm")
# reorder
schaefer400_index_SA$network_label <- factor(schaefer400_index_SA$network_label, levels=c("Vis",
                                                                                          "SomMot", "DorsAttn", "SalVentAttn","Limbic", "Cont", "Default"))
schaefer400_index_SA$network_label[schaefer400_index_SA$network_label=="Limbic"] <- NA
schaefer400_index_SA <- schaefer400_index_SA[order(schaefer400_index_SA$SA12node, 
                                                   schaefer400_index_SA$network_label),]
schaefer400_index_SA$nodeid = c(-1:-376)
tab12 <-table(schaefer400_index_SA$SA12node)
sepid <- cumsum(tab12)
linerange_frame1<-data.frame(y=c(-0.5, -sepid[1:11]+0.5, -sepid[12]-0.5), xmin =rep(0.5, times=13), xmax =rep(1.5, times=13))
linerange_frame2<-data.frame(x=c(0.5,1.5), ymin=c(-376.5, -376.5), ymax=c(-0.5, -0.5))
ggplot()+
  geom_tile(data =schaefer400_index_SA, aes(x=1, y=nodeid, fill = network_label, color=network_label))+
  scale_fill_manual(values=c("#823D8A", "#6E9EC3", "#2E7A3D", "#D867FF", "#DAA34B",
                             "#D76779"))+
  scale_color_manual(values=c("#823D8A", "#6E9EC3", "#2E7A3D", "#D867FF", "#DAA34B",
                              "#D76779"))+
  geom_linerange(data=linerange_frame1, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
  geom_linerange(data=linerange_frame2, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(breaks=NULL, labels = NULL)+
  scale_x_continuous(breaks=NULL, labels = NULL)+
  theme(axis.line = element_blank(), 
        #axis.ticks=element_line(linewidth = 0),
        axis.text.x=element_text(size=12, angle=45, hjust=1), 
        axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
        axis.title =element_text(size=18),
        plot.title = element_text(size=18, hjust = 0.5),
        #legend.position="none",
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(linewidth = 0), 
        panel.grid.minor=element_line(linewidth = 1))
ggsave(paste0(FigureFolder, "/SA12_Yeo7_bin_frequency.tiff"), width=1, height = 12, units = "cm")

# SA12 bin
ggplot()+
  geom_tile(data =schaefer400_index_SA, aes(x=1, y=nodeid, fill = SA12node, color=SA12node))+
  scale_fill_distiller(type="seq", palette = "RdBu",na.value = "grey")+
  scale_color_distiller(type="seq", palette = "RdBu",na.value = "grey")+
  labs(x=NULL, y=NULL)+
  scale_y_continuous(breaks=NULL, labels = NULL)+
  scale_x_continuous(breaks=NULL, labels = NULL)+
  theme(axis.line = element_blank(), 
        #axis.ticks=element_line(linewidth = 0),
        axis.text.x=element_text(size=12, angle=45, hjust=1), 
        axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
        axis.title =element_text(size=18),
        plot.title = element_text(size=18, hjust = 0.5),
        legend.position="none",
        panel.background=element_rect(fill=NA),
        panel.grid.major=element_line(linewidth = 0), 
        panel.grid.minor=element_line(linewidth = 1))
ggsave(paste0(FigureFolder, "/SA12_bin.tiff"), width=1, height = 12, units = "cm")

### schema developmental curve & EF curve
# development
devdf <- readRDS(paste0(interfileFolder, '/plotdata_SA3_forschema.rds'))
devdf2<- devdf[-which(devdf$SClabel %in% c("SCmean.3", "SCmean.4", "SCmean.5")), ]
ggplot(data=devdf2, aes(x=age, y=fit.ratio, group=SClabel, color=SClabel))+
  geom_line(size=1.2, alpha=0.8)+
  scale_color_manual(values=brewer.pal(6, "Set2"))+
  labs(x="Age (8~22)", y="Structural connectivity strength")+
  scale_y_continuous(breaks=NULL, labels = NULL)+
  scale_x_continuous(breaks=NULL, labels = NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=15, color="black"), 
        axis.title =element_text(size=15, color="black"),aspect.ratio = 1,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
ggsave(paste0(FigureFolder, "/SA12_sumSCinvnode_fit/schemadev_curve.svg"), width=12, height = 12, units = "cm")

# cognition
cogdf <- data.frame(age=c(8,22, 8, 22, 8, 22), EF=c(0, 1.2, 0, 1.1, 0, 0.5),
                    SClabel=c("SC1", "SC1", "SC2", "SC2", "SC3", "SC3"))
ggplot(data=cogdf, aes(x=age, y=EF, group=SClabel, color=SClabel))+
  geom_line(size=1.2, alpha=0.8)+
  scale_color_manual(values=brewer.pal(6, "Set2"))+
  labs(x="Age (8~22)", y="Cognition")+
  scale_y_continuous(breaks=NULL, labels = NULL)+
  scale_x_continuous(breaks=NULL, labels = NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=15, color="black"), 
        axis.title =element_text(size=15, color="black"),aspect.ratio = 1,
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
ggsave(paste0(FigureFolder, "/SA12_sumSCinvnode_fit/schemaEF_curve.svg"), width=12, height = 12, units = "cm")




