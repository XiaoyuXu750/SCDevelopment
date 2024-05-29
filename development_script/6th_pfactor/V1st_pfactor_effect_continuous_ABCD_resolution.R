## Sensitivity analysis: different resolutions of large-scale network

library(mgcv)
library(parallel)
library(tidyverse)
library(reshape)
rm(list = ls())
ds.resolution <- 17
elementnum <- ds.resolution*(ds.resolution+1) /2
CVthr=75
wdpath <- getwd()
if (str_detect(wdpath, "Users")){
  resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA', ds.resolution,'/CV', CVthr)
  source(paste0(functionFolder, "/SCrankcorr.R"))
  
}else if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_ABCD'
}

# source function
source(paste0(functionFolder, '/gamminteraction.R'))

detectCores()
# load data
SCdata<-readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatPFactorGeneral.rds'))
SCdata[,c("sex", "handness")] <- lapply(SCdata[,c("sex", "handness")], as.factor)
plotdata <- readRDS(paste0(interfileFolder, '/plotdatasum.df_SA', ds.resolution, '_sumSCinvnode_siteall_CV', CVthr,'.rds'))
# scale the SC strength by their initial strength in the age span
SCdata.diw <- SCdata
for (x in 1:elementnum){
  region <- grep("SC.", names(SCdata), value = T)[x]
  plotdata.tmp <- plotdata[plotdata$SC_label==paste0("SC.", x, "_h"), ]
  SCstrength.diw <- SCdata[,region] / plotdata.tmp$fit[1]
  SCdata.diw[,region] <- SCstrength.diw
}
SCdata.diw[,grep("SC.", names(SCdata), value = T)] <- lapply(SCdata.diw[,grep("SC.", names(SCdata), value = T)], as.numeric)
SCdata.diw[,c("sex", "handness")] <- lapply(SCdata.diw[,c("sex", "handness")], as.factor)

# 1. Interaction effects & p-factor effects
# age by p.factor
dataname <- "SCdata"
smooth_var <- "age" 
int_var.predict.percentile=0.1
covariates <- "sex+mean_fd"
knots=3
set_fx = TRUE
increments = 1000
stats_only=TRUE
int_var <- "GENERAL"
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- mclapply(1:elementnum, function(x){
    region <- grep("SC.", names(SCdata), value=T)[x]
    gamresult <- gamm.smooth.predict.covariateinteraction(region, dataname, smooth_var, int_var, int_var.predict.percentile, covariates, knots, set_fx, increments, stats_only)
    gamresult <- as.data.frame(gamresult)
    
    return(gamresult)
  }, mc.cores = 50)
  gamresult.tmp <- do.call(rbind, resultsum)
  gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
  gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
  gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")
  
  print(paste0(sum(gamresult.tmp$bootstrap_pvalue.fdr<0.05), " edges have significant age by ", int_var, " effect."))
  print(paste0(sum(gamresult.tmp$bootstrap.P.disease.fdr<0.05), " edges have significant ", int_var, " effect."))
  
  saveRDS(gamresult.tmp, paste0(interfileFolder, "/gamresult_Int_age_pFactor_", int_var, "_CV", CVthr,"_SA", ds.resolution,".rds"))
}


# 2. Correlation to S-A connectional axis
SCrank.df.whole <- data.frame(matrix(NA,1,5))
names(SCrank.df.whole) <- c("ds.resolution", "Interest.var",  "r.spearman",  "p.spin", "int_var")
gamresult.tmp <- readRDS(paste0(interfileFolder, "/gamresult_Int_age_pFactor_", int_var, "_CV", CVthr,"_SA", ds.resolution,".rds"))
summary(gamresult.tmp)
gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")
print(paste0(sum(gamresult.tmp$bootstrap_pvalue.fdr<0.05), " edges have significant age by ", int_var, " effect."))
print(paste0(sum(gamresult.tmp$bootstrap.P.disease.fdr<0.05), " edges have significant ", int_var, " effect."))
# SA7: 7 edges have significant general effect.
# SA17: 14 edges have significant general effect.
SCrank.df.age <- SCrankcorr(gamresult.tmp, "IntpartialRsq", ds.resolution, perm.id.full, dsdata=FALSE)
SCrank.df.general <- SCrankcorr(gamresult.tmp, "T.disease", ds.resolution, perm.id.full, dsdata=FALSE)
SCrank.df <- rbind(SCrank.df.age, SCrank.df.general)
SCrank.df
SCrank <- SCrankcorr(gamresult.tmp, "T.disease", ds.resolution, perm.id.full, dsdata=T)$SCrank

# SA7
# ds.resolution  Interest.var r.spearman   p.spearman int_var
# 1             7 IntpartialRsq  0.2980703 0.1234164413 GENERAL
# 2             7     T.disease  0.5955933 0.0008262257 GENERAL

# SA17
# ds.resolution  Interest.var r.spearman   p.spearman
# 1            17 IntpartialRsq  0.2605115 1.144881e-03
# 2            17     T.disease  0.3319914 2.760619e-05

# 3. scatter polts & matrix graph
for (Interest.var in c("IntpartialRsq", "T.disease")){
  tmpvar <- gamresult.tmp[,Interest.var]
  limthr <- max(abs(gamresult.tmp[,Interest.var]), na.rm=T)
  if (Interest.var=="T.disease" & int_var=="general"){
    ytitle=expression(italic("p")*"-factor association ("*italic("T")*" value)")
  }else{ytitle=Interest.var}
  
  # scatter plot
  scatterFig <- ggplot(data=gamresult.tmp)+
    geom_point(aes(x=SCrank, y=tmpvar, color=tmpvar), size=5)+
    geom_smooth(aes(x=SCrank, y=tmpvar),linewidth=1.4, method ="lm", color="black")+
    scale_color_distiller(type="seq", palette = "RdBu", direction = -1, limits=c(-limthr, limthr))+
    labs(x="S-A connectional axis rank", y=ytitle)+
    theme_classic()+
    theme(axis.text=element_text(size=23, color="black"), 
          axis.title =element_text(size=23),aspect.ratio = 0.82,
          axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
          plot.title = element_text(size=20, hjust = 0.5, vjust=2),
          plot.background=element_rect(fill="transparent"),
          panel.background=element_rect(fill="transparent"),
          legend.position = "none")
  scatterFig
  ggsave(paste0(FigureFolder, '/Disease/pFactor/',Interest.var, '_',int_var,'_SCrankcorr.tiff'),scatterFig, width=13, height =13, units = "cm")
  ggsave(paste0(FigureFolder, '/Disease/pFactor/',Interest.var, '_',int_var,'_SCrankcorr.svg'),scatterFig, dpi=600, width=17.5, height =15, units = "cm")
  
  # matrix
  Matrix.tmp.T <- matrix(NA,ds.resolution,ds.resolution)
  Matrix.tmp.T[lower.tri(Matrix.tmp.T, diag = T)] <- tmpvar
  Matrix.tmp.T[upper.tri(Matrix.tmp.T)] <- t(Matrix.tmp.T)[upper.tri(Matrix.tmp.T)]
  colnames(Matrix.tmp.T) <-seq(1, ds.resolution)
  rownames(Matrix.tmp.T) <-seq(1, ds.resolution)
  matrixtmp.df <- as.data.frame(Matrix.tmp.T)
  matrixtmp.df$nodeid <- seq(1, ds.resolution)
  matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
  matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
  matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
  matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
  
  Matrix.tmp.sig <- matrix(NA,ds.resolution,ds.resolution)
  if (Interest.var=="IntpartialRsq"){
    Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (gamresult.tmp$bootstrap_pvalue.fdr < 0.05)
  }else{
    Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (gamresult.tmp$bootstrap.P.disease.fdr < 0.05)
  }
  
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
  titlematrix <- paste0(int_var, "_", Interest.var)
  linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                              y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))
  
  MatFig<-ggplot(data =matrixtmp.df.melt)+
    geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
    scale_fill_distiller(type="seq", palette = "RdBu",na.value = "grey", limits=c(-limthr, limthr))+
    scale_color_distiller(type="seq", palette = "RdBu",na.value = "grey", limits=c(-limthr, limthr))+
    geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "*"), vjust = 0.7, hjust = 0.5, size=6)+
    geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
    geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
    geom_segment(aes(x = 0.5 , y = -0.5 , xend = ds.resolution+0.5 ,yend = -ds.resolution-0.5), color="black", linewidth=0.5)+
    ggtitle(label = titlematrix)+labs(x=NULL, y=NULL)+
    scale_y_continuous(breaks=NULL, labels = NULL)+
    scale_x_continuous(breaks=NULL, labels = NULL)+
    theme(axis.line = element_blank(),
          #axis.ticks=element_line(linewidth = 0),
          axis.text.x=element_text(size=12, angle=45, hjust=1),
          axis.text.y=element_text(size=12, angle=315, hjust=1,vjust=1),
          axis.title =element_text(size=18),
          plot.title = element_text(size=12, hjust = 0.5),
          legend.title=element_text(size=18),
          legend.text=element_text(size=18),
          panel.background=element_rect(fill=NA),
          panel.grid.major=element_line(linewidth = 0),
          panel.grid.minor=element_line(linewidth = 1))
  MatFig
  filename<-paste0(FigureFolder, '/Disease/pFactor/',Interest.var, '_',int_var,'_Matrix', ds.resolution, '.tiff')
  ggsave(filename, MatFig,  height = 18, width = 20, units = "cm")
  
}




