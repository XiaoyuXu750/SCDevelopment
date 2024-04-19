## This script is to fit gam models for each edge in sub-datasets separated at the flip age.
## Statistical indexes and gam model files will be generated.  
library(mgcv)
library(parallel)
rm(list=ls())
CVthr = 75
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
functionFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCP_final/SA12/CV', CVthr)
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.merge.rds'))
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata.diw_CV", CVthr, ".rds"))
summary(SCdata.sum.merge$age)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 8.083  11.500  14.500  14.746  17.875  21.917 
nrow(SCdata.sum.merge)
sepage <- 15.616
SCdata.sum.merge_younger <- SCdata.sum.merge[SCdata.sum.merge$age<sepage,]#n=358 (15.66); n=376 (16)
SCdata.sum.merge_older <- SCdata.sum.merge[SCdata.sum.merge$age>=sepage,]#n=232 (15.66); n=214 (16)
SCdata.younger <- SCdata[SCdata$age<sepage,]
meanlength_younger <- colMeans(SCdata.younger[,which(str_detect(names(SCdata.younger), "length"))])
SCdata.older <- SCdata[SCdata$age>=sepage,]
meanlength_older <- colMeans(SCdata.older[,which(str_detect(names(SCdata.older), "length"))])
meanlength <- colMeans(SCdata[,which(str_detect(names(SCdata), "length"))])

#### source function
source(paste0(functionFolder, '/gamsmooth.R'))
perm.id.full<-readRDS("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA12_sphericalrotations_N10000.rds")
source("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/rotate.parcellation.R")
source("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/perm.sphere.p.R")
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))
detectCores()
## calculate gam results
covariates<-"gender+handnessfactor+race+mean_fd"
smooth_var<-"age"
## younger 8~agesep
dataname<-"SCdata.sum.merge_younger"
resultsum <- mclapply(1:78, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=FALSE, stats_only = FALSE, mod_only=FALSE)
  gamresult<-as.data.frame(gamresult)
  return(gamresult)
}, mc.cores = 3)
gamresultsum.df_younger <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df_younger[,c(2:18)]<-lapply(gamresultsum.df_younger[,c(2:18)], as.numeric)
SCrankcorr(gamresultsum.df_younger, "partialRsq", 12, perm.id.full, dsdata=FALSE)
# CV75
#   ds.resolution matsize Interest.var r.spearman p.spin
# 1            12     376   partialRsq -0.7069528 0.01765
# CV25
# ds.resolution matsize Interest.var r.spearman p.spin
# 1            12     376   partialRsq -0.6932445 0.00055
# control mean length
gamresultsum.df_younger$partialRsq2 <- gamresultsum.df_younger$partialRsq
gamresultsum.df_younger$partialRsq2[which.max(gamresultsum.df_younger$partialRsq2)] <- NA
gamresultsum.df_younger$meanlength <- meanlength_younger
gamresultsum.df_younger$partialRsq_controlmeanlength[!is.na(gamresultsum.df_younger$partialRsq2)] <- residuals(lm(partialRsq2~meanlength, data=gamresultsum.df_younger, na.action = na.omit))
SCrankcorr(gamresultsum.df_younger, "partialRsq_controlmeanlength", 12, perm.id.full, dsdata=FALSE)
#   ds.resolution           Interest.var r.spearman p.spin
# 1            12  partialRsq_controlmeanlength -0.6154746 0.024725

## calculate gam models
resultsum <- mclapply(1:78, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=FALSE, stats_only = TRUE, mod_only=TRUE)
  return(gamresult)
}, mc.cores = 3)

saveRDS(resultsum, paste0(interfileFolder, '/gammodel78_sumSCinvnode_over8_CV',CVthr, '_younger.rds'))
################

## agesep~22
dataname<-"SCdata.sum.merge_older"
resultsum <- mclapply(1:78, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=FALSE, stats_only = FALSE, mod_only=FALSE)
  gamresult<-as.data.frame(gamresult)
  return(gamresult)
}, mc.cores = 3)
gamresultsum.df_older <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
gamresultsum.df_older[,c(2:18)]<-lapply(gamresultsum.df_older[,c(2:18)], as.numeric)
SCrankcorr(gamresultsum.df_older, "partialRsq", 12, perm.id.full, dsdata=FALSE)
# CV75
#   ds.resolution Interest.var r.spearman p.spin
# 1           12 partialRsq    0.4280692 0.013875
# CV25
# ds.resolution matsize Interest.var r.spearman  p.spin
# 1            12     376   partialRsq 0.1214781 0.09675
# when separate age is the same as CV75 (15.66), r=0.22, p=0.0253

# control mean length
gamresultsum.df_older$meanlength <- meanlength_older
gamresultsum.df_older$partialRsq_controlmeanlength <- residuals(lm(partialRsq~meanlength, data=gamresultsum.df_older, na.action = na.omit))
SCrankcorr(gamresultsum.df_older, "partialRsq_controlmeanlength", 12, perm.id.full, dsdata=FALSE)
#   ds.resolution matsize                 Interest.var r.spearman p.spin
# 1            12     376 partialRsq_controlmeanlength 0.363979 0.015225

saveRDS(gamresultsum.df_older, paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV',CVthr, '_older.rds'))

## calculate gam models
resultsum <- mclapply(1:78, function(x){
  SClabel<-names(SCdata.sum.merge)[1+x]
  region<-SClabel
  if (sum(SCdata.sum.merge[,x+1])==0){
    gamresult<-list("NA")
  }else{
    gamresult<-gam.fit.smooth(region, dataname, smooth_var, covariates, knots=3, set_fx=FALSE, stats_only = TRUE, mod_only=TRUE)
  }
  return(gamresult)
}, mc.cores = 3)

saveRDS(resultsum, paste0(interfileFolder, '/gammodel78_sumSCinvnode_over8_CV',CVthr, '_older.rds'))
#######

## younger
### scatter plots
############################################
ds.resolution <- 12
## partial Rsq
gamresultsum.df_younger$partialRsq1 <- gamresultsum.df_younger$partialRsq
gamresultsum.df_younger$partialRsq1[which.max(gamresultsum.df_younger$partialRsq1)] <- NA
computevar <- "partialRsq1"
correlation.df <- SCrankcorr(gamresultsum.df_younger, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
SCrankcorr(gamresultsum.df_younger, computevar, ds.resolution, perm.id.full,dsdata=FALSE)
#   ds.resolution matsize Interest.var r.spearman  p.spin
# 1            12     376  partialRsq1 -0.6832816 0.0067
# CV25 r=-0.6684621 0.00035
summary(gamresultsum.df_younger$partialRsq1)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- gamresultsum.df_younger$partialRsq1
colorbarvalues.Rsq <- colorbarvalues(correlation.df$partialRsq1, 0.2)
RdBucol <- rev(brewer.pal(11, "RdBu"))
RdBu2 <- palette(c(RdBucol[1:6],"#FFFFFF",RdBucol[7], RdBucol[7:9]))
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq1, color=partialRsq1), size=3)+
  geom_smooth(aes(x=SCrank, y=partialRsq1), method ="lm", color="black")+
  #scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.Rsq)+
  scale_colour_gradientn(colours = RdBu2,values=colorbarvalues.Rsq, space="Lab")+
  labs(x="Connectional axis rank", y="Partial R2")+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=15, color="black"), 
        axis.title =element_text(size=15),aspect.ratio = 0.8,
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        axis.line = element_line(linewidth = 0.4),axis.ticks = element_line(linewidth = 0.4),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"), legend.position = "none")
ggsave(paste0(FigureFolder, '/correlation_sumSCinvnode_SCrank_younger/', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/correlation_sumSCinvnode_SCrank_younger/', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=10, height =8, units = "cm")
# matrix plot
computevar <- "partialRsq"
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))
Matrix.tmp <- mtrixplot
Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
colnames(Matrix.tmp) <-seq(1, 12)
rownames(Matrix.tmp) <-seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
colorbarvalues.tmp <- colorbarvalues(matrixtmp.df.melt$value, 0.2)

Matrix.tmp.sig <- matrix(NA, nrow = 12, ncol=12)
Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (gamresultsum.df_younger$pfdr<0.05)
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
ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
  scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = RdBucol[11])+
  scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = RdBucol[11])+
  #geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=9)+
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
filename<-paste0(FigureFolder,"/Matrix12_sumSCinvnode_gamstats_younger/", computevar, "_12net_delLM_CV", CVthr,".tiff")
ggsave(filename, height = 18, width = 20, units = "cm")

## old
### scatter plots
############################################
matsize<-376
ds.resolution <- 12
## partial Rsq
computevar <- "partialRsq"
correlation.df <- SCrankcorr(gamresultsum.df_older, matsize, computevar, ds.resolution, perm.id.full,dsdata=TRUE)
summary(gamresultsum.df_younger$partialRsq)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$partialRsq
colorbarvalues.Rsq <- colorbarvalues(correlation.df$partialRsq, 0.15)
ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=partialRsq, color=partialRsq), size=3)+
  geom_smooth(aes(x=SCrank, y=partialRsq), method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, values=colorbarvalues.Rsq)+
  labs(x="Connectional axis rank", y="Partial R2")+
  #scale_y_continuous(breaks = c(0.0030, 0.0060, 0.0090, 0.012))+
  theme_classic()+
  theme(axis.text=element_text(size=15, color="black"), 
        axis.title =element_text(size=15),aspect.ratio = 0.8,
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        axis.line = element_line(linewidth = 0.4),axis.ticks = element_line(linewidth = 0.4),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"), legend.position = "none")
ggsave(paste0(FigureFolder, '/correlation_sumSCinvnode_SCrank_older/', computevar, '_SCrankcorr_n', ds.resolution, '.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/correlation_sumSCinvnode_SCrank_older/', computevar, '_SCrankcorr_n', ds.resolution, '.svg'), dpi=600, width=10, height =8, units = "cm")
# matrix plot
computevar <- "partialRsq"
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))
Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- correlation.df[,2]
Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
colnames(Matrix.tmp) <-seq(1, 12)
rownames(Matrix.tmp) <-seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
colorbarvalues.tmp <- colorbarvalues(matrixtmp.df.melt$value, 0.15)
ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
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
filename<-paste0(FigureFolder,"/Matrix12_sumSCinvnode_gamstats_older/", computevar, "_12net_delLM_CV", CVthr,".tiff")
ggsave(filename, height = 18, width = 20, units = "cm")


# outputs
#1. gamresults78_over8.rds contains statistic variables from gam models
#2. gammodel78_over8.rds contains 78 gam model files
# younger : 8~15.66 yr
# older : 15.66~22 yr




