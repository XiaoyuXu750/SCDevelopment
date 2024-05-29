# validation for additional covariates.
## This script is to calculate correlation coefficients between cognition var and SC strength
## GAM models were used to regress covariates
library(mgcv)
library(parallel)
library(psych)
library(reshape)
library(RColorBrewer)
library(tidyverse)
rm(list = ls())
CVthr=75
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
addcovariate <- "income.adj" # income.adj / ICV
wdpath <- getwd()
if (str_detect(wdpath, "Users")){
  resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA', ds.resolution,'/CV', CVthr)
  source(paste0(functionFolder, "/SCrankcorr.R"))
  source(paste0(functionFolder, '/colorbarvalue.R'))
}else if (str_detect(wdpath, "cuizaixu_lab")){
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_ABCD'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
}
# load data
SCdata<-readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatNTBfluid_cog2SC.rds'))
# source function
source(paste0(functionFolder, "/gamcog.R"))
detectCores()

## Cognition association
###############
## vars only measured at baseline, gam model was used
Cogvar <- "nihtbx_fluidcomp_uncorrected"
summary(SCdata[,Cogvar])
table(SCdata$eventname[!is.na(SCdata[,Cogvar])])
nonna_index<-which(!is.na(SCdata[ ,Cogvar]))
SCdata.cog<-SCdata[nonna_index,]
SCdata.cog <- SCdata.cog[SCdata.cog$eventname=="baseline_year_1_arm_1",] #n=3836
table(SCdata.cog$siteID)

dataname<-"SCdata.cog"
smooth_var<-"age"
covariates<-paste0("sex+mean_fd+", addcovariate)
knots<-3
corrmethod<-"pearson"

# fit model
summary(SCdata.cog$age)
Nsub<-nrow(SCdata.cog)
if (str_detect(wdpath, "cuizaixu_lab")){
  resultsum <- mclapply(1:elementnum, function(x) {
    SClabel <- grep("SC.", names(SCdata), value=T)[x]
    region <- SClabel
    gamresult <-
      gam.fit.cognition(
        Cogvar,
        dataname,
        region, 
        smooth_var, 
        covariates, 
        knots,
        corrmethod, set_fx = TRUE, stats_only = TRUE
      )
    gamresult <- as.data.frame(gamresult)
    return(gamresult)
  }, mc.cores = 50)
  SC_Cog_results.df <- do.call(rbind, lapply(resultsum, function(x) data.frame(x)))
  SC_Cog_results.df[,c(3:10)]<-lapply(SC_Cog_results.df[,c(3:10)], as.numeric)
  SC_Cog_results.df$corr.p.fdr<-p.adjust(SC_Cog_results.df$corrp, method="fdr")
  SC_Cog_results.df$anova.cov.p.fdr<-p.adjust(SC_Cog_results.df$anova.cov.pvalue, method="fdr")
  SC_Cog_results.df$gam.smooth.p.fdr<-p.adjust(SC_Cog_results.df$gam.smooth.pvalue, method="fdr")
  summary(SC_Cog_results.df)
  SC_Cog_results.df$cognition_var[which.min(SC_Cog_results.df$gam.smooth.t)] # SC.74
  saveRDS(SC_Cog_results.df, paste0(interfileFolder, "/SC_Cog_results_", Cogvar,"_CV", CVthr,"_SA", ds.resolution,"_validat_', addcovariate, '.rds"))
  length(SC_Cog_results.df$cognition_var[SC_Cog_results.df$anova.cov.p.fdr<0.05])}else{
    SC_Cog_results.df <- readRDS(paste0(interfileFolder, "/SC_Cog_results_", Cogvar,"_CV", CVthr,"_SA", ds.resolution,"_validat_', addcovariate, '.rds"))
  }


SC_Cog_results.df.whole <- SC_Cog_results.df
SCrankresult.whole<-SCrankcorr(SC_Cog_results.df.whole,"gam.smooth.t", ds.resolution)
# SES
# ds.resolution Interest.var r.spearman    p.spearman
# 1            12 gam.smooth.t -0.4148161  0.0001593453

# ICV
# ds.resolution Interest.var r.spearman    p.spearman
# 1            12 gam.smooth.t  -0.221774  0.05100872

correlation.df <- SCrankcorr(SC_Cog_results.df, "gam.smooth.t", 12,dsdata=TRUE)
correlation.df$sig <- (SC_Cog_results.df$anova.cov.p.fdr<0.05)
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))
SC_Cog_results.tmp <-SC_Cog_results.df.whole
SC_Cog_results.tmp$SCrank <- correlation.df$SCrank

lwth <- min(SC_Cog_results.tmp$gam.smooth.t)
Fig <- ggplot(data=SC_Cog_results.tmp)+
  geom_point(aes(x=SCrank, y=gam.smooth.t, color=gam.smooth.t), size=5)+
  geom_smooth(aes(x=SCrank, y=gam.smooth.t), method ="lm", color="black", linewidth=1.4)+
  scale_colour_distiller(type="seq", palette = "RdBu",limits=c(lwth, -lwth), direction = -1)+
  labs(x="S-A connectional axis rank", y=expression("Cognitive association ("*italic("T")*" value)"))+
  theme_classic()+theme_classic()+theme(axis.text=element_text(size=23, color="black"), 
                                        axis.title =element_text(size=23),aspect.ratio = 0.82,
                                        axis.line = element_line(linewidth=0.6),axis.ticks= element_line(linewidth=0.6),
                                        plot.title = element_text(size=20, hjust = 0.5, vjust=2),
                                        plot.background=element_rect(fill="transparent"),
                                        panel.background=element_rect(fill="transparent"),
                                        legend.position = "none")

Fig
ggsave(paste0(FigureFolder, '/cognition/', Cogvar, '/CorrTvalue_SCrankcorr_n12_siteall_', addcovariate, '.tiff'), Fig, width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/cognition/', Cogvar, '/CorrTvalue_SCrankcorr_n12_siteall_', addcovariate, '.svg'), Fig, dpi=600, width=17.5, height =15, units = "cm")






