# validation for different matrix resolution: 17 and 7
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
ds.resolution <- 7
elementnum <- ds.resolution*(ds.resolution+1) /2
wdpath <- getwd()
if (str_detect(wdpath, "Users")){
  resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA', ds.resolution,'/CV', CVthr)
  source(paste0(functionFolder, "/SCrankcorr.R"))
  perm.id.full<-readRDS(paste0("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA", ds.resolution,"_sphericalrotations_N10000.rds"))
  source(paste0(functionFolder, '/colorbarvalue.R'))
}else if (str_detect(wdpath, "cuizaixu_lab")){
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_ABCD'
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
}
# load data
SCdata<-readRDS(paste0(interfileFolder, '/SCdata_SA', ds.resolution, '_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatNTBfluid.rds'))
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

dataname<-"SCdata.cog" #n=3836, Male 2019, age 8.92~11.00
smooth_var<-"age"
covariates<-"gender+handness+race+mean_fd"
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
  saveRDS(SC_Cog_results.df, paste0(interfileFolder, "/SC_Cog_results_", Cogvar,"_CV", CVthr,"_SA", ds.resolution,".rds"))
  length(SC_Cog_results.df$cognition_var[SC_Cog_results.df$anova.cov.p.fdr<0.05])}else{
    SC_Cog_results.df <- readRDS(paste0(interfileFolder, "/SC_Cog_results_", Cogvar,"_CV", CVthr,"_SA", ds.resolution,".rds"))
  }

# fluid, 60 edges sig after fdr (SA17); 21 (SA7)

SC_Cog_results.df.whole <- SC_Cog_results.df
SCrankresult.whole<-SCrankcorr(SC_Cog_results.df.whole,"gam.smooth.t", ds.resolution, perm.id.full)
# fluid, r=-0.4360985, p=0.0239, site all (SA17); r=-0.7064459 p=0.0111, (SA7)

### plot: scatter plot+matrix graph
########################
correlation.df <- SCrankcorr(SC_Cog_results.df, "gam.smooth.t", ds.resolution, perm.id.full,dsdata=TRUE)
correlation.df$sig <- (SC_Cog_results.df$anova.cov.p.fdr<0.05)
Matrix.tmp <- matrix(NA, nrow = ds.resolution, ncol=ds.resolution)
linerange_frame<-data.frame(x=c(0.5,ds.resolution+0.5), ymin =rep(-ds.resolution-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -ds.resolution-0.5), xmin=rep(0.5, times=2), xmax=rep(ds.resolution+0.5, times=2))
SC_Cog_results.tmp <-SC_Cog_results.df.whole
SC_Cog_results.tmp$SCrank <- correlation.df$SCrank

lwth <- min(SC_Cog_results.tmp$gam.smooth.t) # 5.394 (SA17); 5.299 (SA7)
Fig <- ggplot(data=SC_Cog_results.tmp)+
  geom_point(aes(x=SCrank, y=gam.smooth.t, color=gam.smooth.t), size=5)+
  geom_smooth(aes(x=SCrank, y=gam.smooth.t), method ="lm", color="black", linewidth=1.4)+
  scale_colour_distiller(type="seq", palette = "RdBu",limits=c(lwth, -lwth), direction = -1)+
  labs(x="Connectional axis rank", y=expression("Cognitive association ("*italic("T")*" value)"))+
  theme_classic()+
  theme(axis.text=element_text(size=23.2, color="black"), 
        axis.title =element_text(size=23.2),aspect.ratio = 1.03,axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(linewidth = 0.6),
        plot.title = element_text(size=15, hjust = 0.5, vjust=0),
        plot.subtitle = element_text(size=15, hjust = 0.9, vjust=-6),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),plot.margin = margin(t=10,r=5,b=5,l=5, unit="pt"),
        legend.position = "none")
Fig
ggsave(paste0(FigureFolder, '/cognition/', Cogvar, '_rev/CorrTvalue_SCrankcorr_n', ds.resolution, '_siteall.tiff'), Fig, width=17, height =14, units = "cm")
ggsave(paste0(FigureFolder, '/cognition/', Cogvar, '_rev/CorrTvalue_SCrankcorr_n', ds.resolution, '_siteall.svg'), Fig, dpi=600, width=14.8, height =13.8, units = "cm")

Matrix.tmp.T <- Matrix.tmp
Matrix.tmp.T[lower.tri(Matrix.tmp.T, diag = T)] <- SC_Cog_results.tmp$gam.smooth.t
Matrix.tmp.T[upper.tri(Matrix.tmp.T)] <- t(Matrix.tmp.T)[upper.tri(Matrix.tmp.T)]
colnames(Matrix.tmp.T) <-seq(1, ds.resolution)
rownames(Matrix.tmp.T) <-seq(1, ds.resolution)
matrixtmp.df <- as.data.frame(Matrix.tmp.T)
matrixtmp.df$nodeid <- seq(1, ds.resolution)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)

Matrix.tmp.sig <- Matrix.tmp
Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- correlation.df$sig
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
titlematrix <- Cogvar
Fig<-ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
  scale_fill_distiller(type="seq", palette = "RdBu",limits=c(lwth, -lwth),na.value = "grey")+
  scale_color_distiller(type="seq", palette = "RdBu",limits=c(lwth, -lwth),na.value = "grey")+
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
Fig
filename<-paste0(FigureFolder, '/cognition/', Cogvar, '_rev/CorrTvalue_Matrix_n', ds.resolution, '_siteall.tiff')
ggsave(filename, Fig,  height = 18, width = 20, units = "cm")




