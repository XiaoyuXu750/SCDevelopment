library(mgcv)
library(psych)
library(tidyverse)
library(parallel)
library(scales)
library(openxlsx)
library(gratia)
library(ggplot2)
library(RColorBrewer)
library(paletteer)
rm(list = ls())
# Set path
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results_HCPD'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'
# Load data
CVthr = 75
gamresult <- readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatage.rds'))
derivative <- readRDS(paste0(resultFolder, '/derivative.df78_CV', CVthr,'.rds'))
plotdatasum <- readRDS(paste0(interfileFolder, '/plotdatasum_scale_TRUE_SA12.rds'))
SA12_otherorg <- read.csv(paste0(interfileFolder, '/SA12_OrtherOrganization.csv'))

gammod.younger <- readRDS(paste0(interfileFolder, '/gammodel78_sumSCinvnode_over8_CV',CVthr, '_younger.rds'))
gammod.older <- readRDS(paste0(interfileFolder, '/gammodel78_sumSCinvnode_over8_CV',CVthr, '_older.rds'))


source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/plotdata_generate.R'))
source('D:/xuxiaoyu/DMRI_network_development/Normative_model/functions/plotmatrix.R')

### 1. Competing axis
######################
matrix.G1 <- matrix(NA, 12, 12)
matrix.T1T2ratio <- matrix(NA, 12, 12)
for (x in 1:12){
  for (y in 1:12){
    matrix.G1[x, y] <- (SA12_otherorg$MeanG1[x]/abs(SA12_otherorg$MeanG1[x]))*SA12_otherorg$MeanG1[x] ^2 + 
      (SA12_otherorg$MeanG1[y]/abs(SA12_otherorg$MeanG1[y]))*SA12_otherorg$MeanG1[y] ^2
    
    matrix.T1T2ratio[x, y] <- (SA12_otherorg$MeanT1T2ratio[x]/abs(SA12_otherorg$MeanT1T2ratio[x]))*SA12_otherorg$MeanT1T2ratio[x] ^2 + 
      (SA12_otherorg$MeanT1T2ratio[y]/abs(SA12_otherorg$MeanT1T2ratio[y]))*SA12_otherorg$MeanT1T2ratio[y] ^2
  }
}
indexlowtri <- lower.tri(matrix.G1, diag = T)
OrtherAxis <- data.frame(SClabel=paste0("SC.", 1:78), G1 = rank(matrix.G1[indexlowtri]), 
                         T1T2ratio = rank(matrix.T1T2ratio[indexlowtri]))
PaletteSet <- list(Name="RdBu", drirection=-1, lmmin = 1, lmmax = 78, 
                   anglex=90, angley=0, hjustx=1, vjustx=0.5, hjusty=1, vjusty=0.5)
Mat12.G1 <- plotmatrix(dataname = "OrtherAxis", variable="G1", ds.resolution=12, Pvar=NA, NAcol="white", lmthr=NA, 
                       axeslabels=NULL, axeslabelsGap=T, linerange_frame=NA, PaletteSet=PaletteSet)
ggsave(paste0(FigureFolder, '/CV', CVthr, '/PCA_corr/G1fMRI_axis.tiff'), Mat12.G1,  height = 18, width = 20, units = "cm")


Mat12.T1T2ratio <- plotmatrix(dataname = "OrtherAxis", variable="T1T2ratio", ds.resolution=12, Pvar=NA, NAcol="white", lmthr=NA, 
                              axeslabels=NULL, axeslabelsGap=T, linerange_frame=NA, PaletteSet=PaletteSet)
ggsave(paste0(FigureFolder, '/CV', CVthr, '/PCA_corr/T1T2ratio_axis.tiff'), Mat12.T1T2ratio,  height = 18, width = 20, units = "cm")

SCrankcorr(OrtherAxis, "G1", 12) # rho = 0.9846224, P=2.3999e-59
SCrankcorr(OrtherAxis, "T1T2ratio", 12) # rho = -0.9043958 P=7.787875e-30

######################


# PCA on 1st derivative
######################################
derivative <- as.data.frame(derivative)
derivative.df <- derivative %>% select(label_ID, derivative, age) %>% pivot_wider(names_from = "label_ID", values_from = "derivative")

derivative.df$age <- NULL

derivative.pca <- prcomp(derivative.df, rank. = 3)
summary(derivative.pca)
# Importance of first k=3 (out of 78) components:
#   PC1       PC2       PC3
# Standard deviation     0.05489 3.543e-10 4.658e-17
# Proportion of Variance 1.0000 0.000e+00 0.000e+00
# Cumulative Proportion  1.0000 1.000e+00 1.000e+00


## All
derivative.df$age <- NULL
derivative.PCA.df <- data.frame(SClabel = names(derivative.df), PC1.loading =derivative.pca$rotation[,1])
derivative.PCA.df$index <- gsub("SC.", "", derivative.PCA.df$SClabel)
derivative.PCA.df$index <- gsub("_h", "", derivative.PCA.df$index)
derivative.PCA.df <- derivative.PCA.df[order(as.numeric(derivative.PCA.df$index)), ]
SCrankcorr(derivative.PCA.df, "PC1.loading", 12)
# ds.resolution Interest.var r.spearman p.spearman
# 1            12  PC1.loading  0.8030882 9.189684e-19

derivative.PCA.df$SCrank <- SCrankcorr(derivative.PCA.df, "PC1.loading", 12, T)$SCrank

resultdf$SAaxis.corr[2] <- SCrankcorr(derivative.PCA.df, "PC1.loading", 12)$r.spearman
resultdf$SAaxis.P[2] <- SCrankcorr(derivative.PCA.df, "PC1.loading", 12)$p.spearman

## Correlation with other axis
resultdf$G1axis.corr[2] <- corr.test(derivative.PCA.df$PC1.loading, OrtherAxis$G1, method = "spearman")$r
resultdf$G1axis.P[2] <- corr.test(derivative.PCA.df$PC1.loading, OrtherAxis$G1, method = "spearman")$p

resultdf$T1T2axis.corr[2] <- corr.test(derivative.PCA.df$PC1.loading, OrtherAxis$T1T2ratio, method = "spearman")$r
resultdf$T1T2axis.P[2] <- corr.test(derivative.PCA.df$PC1.loading, OrtherAxis$T1T2ratio, method = "spearman")$p

# plot
# S-A axis
ggplot(data=derivative.PCA.df)+
  geom_point(aes(x=SCrank, y=PC1.loading, color=SCrank), size=5)+
  geom_smooth(aes(x=SCrank, y=PC1.loading),linewidth=2, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="S-A connectional axis rank", y="PC1 loadings of \nchange-rate trajectories")+
  theme_classic()+mytheme

ggsave(paste0(FigureFolder,'/CV', CVthr, '/PCA_corr/PC1_SCrank_1stderivTrajectory.tiff'), width=width, height=height, units = "cm")


derivative.PC1 <- plotmatrix(dataname = "derivative.PCA.df", variable="PC1.loading", ds.resolution=12, Pvar=NA, NAcol="white", lmthr=NA, 
                      axeslabels=NULL, axeslabelsGap=T, linerange_frame=NA, PaletteSet=NA)
derivative.PC1
filename<-paste0(FigureFolder, '/CV', CVthr, '/PCA_corr/Mat_derivative_PC1_loading.tiff')
ggsave(filename, derivative.PC1,  height = 18, width = 20, units = "cm")


# G1
derivative.PCA.df$G1_Axis <- OrtherAxis$G1
ggplot(data=derivative.PCA.df)+
  geom_point(aes(x=G1_Axis, y=PC1.loading, color=G1_Axis), size=5)+
  geom_smooth(aes(x=G1_Axis, y=PC1.loading),linewidth=2, method ="lm", color="black")+
  paletteer::scale_color_paletteer_c("pals::coolwarm")+
  labs(x="Gradient1 connectional axis rank", y="PC1 loadings of \nchange-rate trajectories")+
  theme_classic()+mytheme

ggsave(paste0(FigureFolder,'/CV', CVthr, '/PCA_corr/PC1_FCG1_1stderivTrajectory.tiff'), width=width, height=height, units = "cm")

# T1T2ratio
derivative.PCA.df$T1T2ratio_Axis <- OrtherAxis$T1T2ratio
ggplot(data=derivative.PCA.df)+
  geom_point(aes(x=T1T2ratio_Axis, y=PC1.loading, color=T1T2ratio_Axis), size=5)+
  geom_smooth(aes(x=T1T2ratio_Axis, y=PC1.loading),linewidth=2, method ="lm", color="black")+
  paletteer::scale_color_paletteer_c("pals::coolwarm")+
  labs(x="T1/T2 ratio connectional axis rank", y="PC1 loadings of \nchange-rate trajectories")+
  theme_classic()+mytheme

ggsave(paste0(FigureFolder,'/CV', CVthr, '/PCA_corr/PC1_T1T2ratio_1stderivTrajectory.tiff'), width=width, height=height, units = "cm")



