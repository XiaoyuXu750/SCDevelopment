library(mgcv)
library(parallel)
library(tidyverse)
library(reshape)
rm(list = ls())
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
CVthr=75
wdpath <- getwd()
addcovariate <- "income.adj" # income.adj / ICV
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
covariates <- paste0("sex+mean_fd+", addcovariate)
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
  
  saveRDS(gamresult.tmp, paste0(interfileFolder, "/gamresult_Int_age_pFactor_", int_var, "_CV", CVthr,"_SA", ds.resolution,"_validat_", addcovariate, ".rds"))
}


# 2. Correlation to S-A connectional axis
SCrank.df.whole <- data.frame(matrix(NA,1,5))
names(SCrank.df.whole) <- c("ds.resolution", "Interest.var",  "r.spearman",  "p.spin", "int_var")
gamresult.tmp <- readRDS(paste0(interfileFolder, "/gamresult_Int_age_pFactor_", int_var, "_CV", CVthr,"_SA", ds.resolution,"_validat_", addcovariate, ".rds"))
summary(gamresult.tmp)
gamresult.tmp[3:12] <- lapply(gamresult.tmp[3:12], as.numeric)
gamresult.tmp$bootstrap_pvalue.fdr <- p.adjust(gamresult.tmp$bootstrap_pvalue, method = "fdr")
gamresult.tmp$bootstrap.P.disease.fdr <- p.adjust(gamresult.tmp$bootstrap.P.disease, method = "fdr")
print(paste0(sum(gamresult.tmp$bootstrap_pvalue.fdr<0.05), " edges have significant age by ", int_var, " effect."))
print(paste0(sum(gamresult.tmp$bootstrap.P.disease.fdr<0.05), " edges have significant ", int_var, " effect."))

SCrank.df.age <- SCrankcorr(gamresult.tmp, "IntpartialRsq", ds.resolution, perm.id.full, dsdata=FALSE)
SCrank.df.general <- SCrankcorr(gamresult.tmp, "T.disease", ds.resolution, perm.id.full, dsdata=FALSE)
SCrank.df <- rbind(SCrank.df.age, SCrank.df.general)
SCrank.df$int_var <- int_var
SCrank.df

# SES
# ds.resolution  Interest.var r.spearman    p.spearman
# 1            12 IntpartialRsq  0.1100081  3.376629e-01
# 2            12     T.disease  0.4820426  7.872256e-06

# ICV
# ds.resolution  Interest.var r.spearman     p.spearman
# 1            12 IntpartialRsq  0.2630507   1.997415e-02
# 2            12     T.disease  0.5013911   2.926239e-06
limthr <- max(abs(gamresult.tmp$T.disease), na.rm=T)
ytitle=expression(italic("p")*"-factor association ("*italic("T")*" value)")
SCrank.tmp <- SCrankcorr(gamresult.tmp, "T.disease", 12, perm.id.full, dsdata=T)
SCrank <- SCrank.tmp$SCrank
scatterFig <- ggplot(data=gamresult.tmp)+
  geom_point(aes(x=SCrank, y=T.disease, color=T.disease), size=5)+
  geom_smooth(aes(x=SCrank, y=T.disease),linewidth=1.4, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, limits=c(-limthr, limthr))+
  labs(x="S-A connectional axis rank", y=ytitle)+
  theme_classic()+
  theme(axis.text=element_text(size=22.6, color="black"), 
        axis.title =element_text(size=22.6),aspect.ratio = 0.82,
        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
scatterFig
ggsave(paste0(FigureFolder, '/Disease/pFactor/T.disease_',int_var,'_SCrankcorr_', addcovariate, '.tiff'),scatterFig, width=13, height =13, units = "cm")
ggsave(paste0(FigureFolder, '/Disease/pFactor/T.disease_',int_var,'_SCrankcorr_', addcovariate, '.svg'),scatterFig, dpi=600,width=17.5, height =15, units = "cm")
