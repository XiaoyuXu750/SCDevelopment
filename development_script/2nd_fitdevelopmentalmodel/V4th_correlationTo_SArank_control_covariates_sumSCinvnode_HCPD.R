## Validation: covariates SES / ICV
#### This script is to conduct correlation analysis 
#### between gam statistical indexes to S-A connectional axis rank.
#### And draw scatter plots & matrix graphs. Fig. S10 & S11.
library(R.matlab)
library(tidyverse)
library(parallel)
library(psych)
library(corrplot)
library(reshape)
rm(list = ls())
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
addcovariate <- "income.adj" # income.adj / ICV
demopath<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/demopath'
functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
resultFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results'
interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final/SA', ds.resolution)

#### load data
CVthr = 75
gamresult<-readRDS(paste0(interfileFolder, '/gamresults', elementnum, '_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE_validat_',addcovariate,'.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.merge.rds'))
#### source function
source(paste0(functionFolder, '/SCrankcorr.R'))
source(paste0(functionFolder, '/colorbarvalue.R'))

#### convert critical ages of insignificantly developmental edges to NA
#### convert critical ages equal to age boundaries to NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method="fdr")
gamresult$sig <- (gamresult$pfdr < 0.05)
gamresult$increase.onset[gamresult$sig==FALSE]<-NA ; gamresult$increase.onset2 <- gamresult$increase.onset
gamresult$increase.onset2[round(gamresult$increase.onset2,2)==8.08] <- NA
gamresult$increase.offset[gamresult$sig==FALSE]<-NA ; gamresult$increase.offset2 <- gamresult$increase.offset
gamresult$increase.offset2[round(gamresult$increase.offset2,2)==21.92] <- NA
gamresult$peak.change[gamresult$sig==FALSE]<-NA
gamresult$peak.increase.change[gamresult$sig==FALSE]<-NA
gamresult$pfdr <- p.adjust(gamresult$anova.smooth.pvalue, method = "fdr")

#### 1. compute correlations to SC rank
computevar.list <- c("partialRsq", "increase.onset2", "increase.offset2", "peak.increase.change",
                     "meanderv2")
for (x in 1:5){
  computevar <- computevar.list[x]
  ds.resolution<-ds.resolution
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=FALSE)
  if (x==1){
    SCrank_correlation <- correlation.df
  }else{
    SCrank_correlation <- rbind(SCrank_correlation, correlation.df)
  }
}
SCrank_correlation
# total strength
# ds.resolution         Interest.var r.spearman   p.spin   p.spearman
# 1            12           partialRsq -0.3857934 0.103575 4.860970e-04
# 2            12      increase.onset2  0.7298701 0.002675 1.729791e-04
# 3            12     increase.offset2  0.4965035 0.076750 1.006026e-01
# 4            12 peak.increase.change  0.7755106 0.004150 1.252191e-06
# 5            12            meanderv2  0.7932242 0.002475 4.837113e-18

# ICV
# ds.resolution         Interest.var r.spearman   p.spin   p.spearman
# 1            12           partialRsq -0.3599069 0.100925 1.210138e-03
# 2            12      increase.onset2  0.1704853 0.233300 3.765852e-01
# 3            12     increase.offset2  0.4156952 0.023500 1.441821e-03
# 4            12 peak.increase.change  0.6999653 0.000950 2.192115e-11
# 5            12            meanderv2  0.8405332 0.002600 6.347451e-22

# income.adj
# ds.resolution         Interest.var r.spearman   p.spin   p.spearman
# 1            12           partialRsq -0.2164879 0.167600 5.694358e-02
# 2            12      increase.onset2  0.3359578 0.025725 2.405314e-02
# 3            12     increase.offset2  0.6001861 0.004075 6.800919e-05
# 4            12 peak.increase.change  0.8242220 0.000250 3.324529e-17
# 5            12            meanderv2  0.7951338 0.002875 3.532107e-18

### 2. scatter plots
############################################
ds.resolution <- ds.resolution
# meanderv2
computevar <- "meanderv2"
correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
mtrixplot <- matrix(NA, ds.resolution, ds.resolution)
mtrixplot[lower.tri(mtrixplot, diag = T)] <- correlation.df$meanderv2

ggplot(data=correlation.df)+
  geom_point(aes(x=SCrank, y=meanderv2, color=SCrank), size=5.5)+
  geom_smooth(aes(x=SCrank, y=meanderv2),linewidth=2, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="S-A connectional axis rank", y="Second derivative")+
  scale_y_continuous(breaks = c(-0.002, 0, 0.002), labels = c(-2, 0, 2))+
  theme_classic()+
  theme(axis.text=element_text(size=23, color="black"), 
        axis.title =element_text(size=23),aspect.ratio = 0.75,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '_',addcovariate,'.tiff'), width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV', CVthr, '/correlation_sumSCinvnode_SCrank/mean', computevar, '_SCrankcorr_n', ds.resolution, '_',addcovariate,'.svg'), dpi=600, width=17.5, height =15, units = "cm")

### 3. matrix graphs for resolution of ds.resolution
############################################
Matrix.tmp <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
computevar.list <- c("partialRsq", "increase.onset", "increase.offset", "peak.increase.change", "meanderv2")
linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))
SCrank_correlation.df<-mclapply(1:5, function(x){
  computevar <- computevar.list[x]
  ds.resolution<-ds.resolution
  correlation.df <- SCrankcorr(gamresult, computevar, ds.resolution, dsdata=TRUE)
  return(correlation.df)
}, mc.cores = 4)
colorbar.prob <- c(0.5, 0.4, 0.6, 0.5, 0.44, 0.6)

for (i in 1:5){
  computevar <- computevar.list[i]
  Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank_correlation.df[[i]][,2]
  Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
  colnames(Matrix.tmp) <-seq(1, ds.resolution)
  rownames(Matrix.tmp) <-seq(1, ds.resolution)
  matrixtmp.df <- as.data.frame(Matrix.tmp)
  matrixtmp.df$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  colorbarvalues.tmp <- colorbarvalues(SCrank_correlation.df[[i]][,2], colorbar.prob[i])
  if (computevar=="partialRsq"){
    lmthr <- max(abs(gamresult$partialRsq))
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
    
    Fig<-ggplot(data =matrixtmp.df.melt)+
      geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
      scale_fill_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
      scale_color_distiller(type="seq", palette = "RdBu",limit=c(-lmthr, lmthr), na.value = "grey")+
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
      #paletteer::scale_fill_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.tmp,oob = squish) +
      #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1,values=colorbarvalues.tmp, oob = squish) +
      scale_fill_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "grey")+
      scale_color_distiller(type="seq", palette = "RdBu",values=colorbarvalues.tmp, na.value = "grey")+
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
  filename<-paste0(FigureFolder,'/CV', CVthr,  "/Matrix", ds.resolution, "_sumSCinvnode_gamstats_Age8_22/", computevar, "_", ds.resolution, "net_delLM_CV75_",addcovariate,".tiff")
  ggsave(filename, Fig,  height = 18, width = 20, units = "cm")
  
}



