## This script is to fit interaction models between ADHD and HC.
## Both age and cognition are defined as the variable of interest.
## Gamm models were used to regress covariables.
## Comorbidity is not considered.
library(mgcv)
library(parallel)
library(tidyverse)
library(reshape)
rm(list = ls())
CVthr=75
wdpath <- getwd()
if (str_detect(wdpath, "Users")){
  resultFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/results_ABCD'
  interfileFolder <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-paste0('/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final/SA12/CV', CVthr)
  source(paste0(functionFolder, "/SCrankcorr_beforegam.R"))
  perm.id.full<-readRDS("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA12_sphericalrotations_N10000.rds")
  
}else if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/Rcode_SCdevelopment/gamfunction'
  resultFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SCdevelopment/results_ABCD'
}

# source function
source(paste0(functionFolder, '/linearmediation.R'))

# load data
SCdata<-readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatPFactorGeneral.rds'))
SCdata[,c("gender", "handness", "race")] <- lapply(SCdata[,c("gender", "handness", "race")], as.factor)
Cogdf <- SCdata %>% select(subID, eventname, nihtbx_fluidcomp_uncorrected) %>% drop_na() %>% 
  filter(str_detect(eventname, "base"))
Cogdf <- reshape::rename(Cogdf, c(nihtbx_fluidcomp_uncorrected="nihtbx_fluidcomp_uncorrected_base"))
Cogdf$eventname <- NULL
SCdata <- SCdata %>% left_join(Cogdf, by="subID")

SCdata.base <- SCdata %>% filter(eventname=="baseline_year_1_arm_1") %>% drop_na(nihtbx_fluidcomp_uncorrected_base)
SCdata.2year <- SCdata %>% filter(eventname=="2_year_follow_up_y_arm_1") %>% drop_na(nihtbx_fluidcomp_uncorrected_base)
interest_vars <- c("general", "nihtbx_fluidcomp_uncorrected_base")
corr.test(SCdata.base[,interest_vars]) # r=-0.14, p=0
corr.test(SCdata.2year[,interest_vars]) # r=-0.1, p=0

# 1. mediation model, cross-section, baseline
# p-factor --> SC --> cognition
dataname <- "SCdata.base"
Y.var <- "nihtbx_fluidcomp_uncorrected_base"
X.var <- "general"
covariates <- "s(age, k=3, fx=T)+gender+handness+race+mean_fd"

# resultsum <- mclapply(1:78, function(x){
#   M.var <- grep("SC.", names(SCdata), value=T)[x]
#   med.result <- linearmediation(dataname, M.var, Y.var, X.var, covariates)
#   med.result <- as.data.frame(med.result)
# 
#   return(med.result)
# }, mc.cores=50)
# saveRDS(resultsum, paste0(resultFolder, '/mediation_', X.var, '_SC_', Y.var, '_CV', CVthr, '.rds'))
resultsum <- readRDS(paste0(resultFolder, '/mediation_', X.var, '_SC_', Y.var, '_CV', CVthr, '.rds'))
med.base <- do.call(rbind, resultsum)
med.base[,3:22] <- lapply(med.base[,3:22], as.numeric)
med.base$indirect.P.fdr <- p.adjust(med.base$indirect.P, method = "fdr")
summary(med.base$indirect.P.fdr) # 0 edge sig
med.base$M.var[med.base$indirect.P < 0.05] # 16 edges sig, before fdr.

SCrank.med <- SCrankcorr(med.base, "indirect.Z", 12, perm.id.full, dsdata=FALSE)
#   ds.resolution Interest.var r.spearman p.spin
# 1            12   indirect.Z -0.4660327 0.0515
SCrank.med.df <- SCrankcorr(med.base, "indirect.Z", 12, perm.id.full, dsdata=T)

# scatter plot
limthr <- max(abs(SCrank.med.df$indirect.Z)) #2.89956
scatterFig <- ggplot(data=SCrank.med.df)+
  geom_point(aes(x=SCrank, y=indirect.Z, color=indirect.Z), size=4)+
  geom_smooth(aes(x=SCrank, y=indirect.Z),linewidth=1.4, method ="lm", color="black")+
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1, limits=c(-limthr, limthr))+
  labs(x="Connectional axis rank", y=expression("Indirect effect ("*italic("Z")*")"))+
  theme_classic()+
  theme(axis.text=element_text(size=20.2, color="black"), 
        axis.title =element_text(size=20.2),aspect.ratio = 1,
        plot.title = element_text(size=21, hjust = 0.5, vjust=2),
        plot.background=element_rect(fill="transparent"),
        panel.background=element_rect(fill="transparent"),
        legend.position = "none")
scatterFig
ggsave(paste0(FigureFolder, '/Disease/pFactor/general_SC_NIHfluid_med_SCrankcorr.tiff'),scatterFig, width=13, height =13, units = "cm")
ggsave(paste0(FigureFolder, '/Disease/pFactor/general_SC_NIHfluid_med_SCrankcorr.svg'),scatterFig, width=13, height =13, units = "cm")

# matrix
Matrix.tmp.T <- matrix(NA,12,12)
Matrix.tmp.T[lower.tri(Matrix.tmp.T, diag = T)] <- SCrank.med.df$indirect.Z
Matrix.tmp.T[upper.tri(Matrix.tmp.T)] <- t(Matrix.tmp.T)[upper.tri(Matrix.tmp.T)]
colnames(Matrix.tmp.T) <-seq(1, 12)
rownames(Matrix.tmp.T) <-seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp.T)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)

Matrix.tmp.sig <- matrix(NA,12,12)
Matrix.tmp.sig[lower.tri(Matrix.tmp.sig, diag = T)] <- (med.base$indirect.P < 0.05)

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
titlematrix <- expression("Indirect effect ("*italic("Z")*")")
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))

MatFig<-ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
  scale_fill_distiller(type="seq", palette = "RdBu",na.value = "grey", limits=c(-limthr, limthr))+
  scale_color_distiller(type="seq", palette = "RdBu",na.value = "grey", limits=c(-limthr, limthr))+
  #scale_colour_gradientn(colours = RdBucol,limits=c(lwth, -lwth), space="Lab")+
  #scale_fill_gradientn(colours = RdBucol,limits=c(lwth, -lwth), space="Lab")+
  geom_text(data =matrixtmp.df.sig.melt, aes(x=variable, y=nodeid, label = "."), vjust = 0.1, hjust = 0.5, size=10)+
  geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
  geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
  geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
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
filename<-paste0(FigureFolder, '/Disease/pFactor/general_SC_NTBfluid_med_IndirectZ_Matrix12.tiff')
ggsave(filename, MatFig,  height = 18, width = 20, units = "cm")


# 2. mediation model, longitudinal, baseline -- 2year
# cognition --> SC --> p-factor
dataname <- "SCdata.2year"
Y.var <- "general"
X.var <- "nihtbx_fluidcomp_uncorrected_base"
covariates <- "s(age, k=3, fx=T)+gender+handness+race+mean_fd"

# resultsum <- mclapply(1:78, function(x){
#   M.var <- grep("SC.", names(SCdata), value=T)[x]
#   med.result <- linearmediation(dataname, M.var, Y.var, X.var, covariates)
#   med.result <- as.data.frame(med.result)
#   
#   return(med.result)
# }, mc.cores=50)
# saveRDS(resultsum, paste0(resultFolder, '/mediation_base', X.var, '_2YSC_2Y', Y.var, '_CV', CVthr, '.rds'))
resultsum <- readRDS(paste0(resultFolder, '/mediation_base', X.var, '_2YSC_2Y', Y.var, '_CV', CVthr, '.rds'))
med.2year <- do.call(rbind, resultsum)
med.2year[,3:22] <- lapply(med.2year[,3:22], as.numeric)
med.2year$indirect.P.fdr <- p.adjust(med.2year$indirect.P, method = "fdr")
summary(med.2year$indirect.P.fdr) # 0 edge sig
med.2year$M.var[med.2year$indirect.P < 0.05] # 5 edges sig, before fdr.

SCrank.med <- SCrankcorr(med.2year, "indirect.Z", 12, perm.id.full, dsdata=FALSE)
#   ds.resolution Interest.var r.spearman p.spin
# 1            12   indirect.Z -0.2827407 0.134175

















