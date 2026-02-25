# Test within-person effects
rm(list=ls())
library(mgcv)
library(parallel)
library(tidyverse)
library(lme4)
library(lmerTest)
wdpath <- getwd()
# set resolution
ds.resolution <- 12
elementnum <- ds.resolution*(ds.resolution+1) /2
if (str_detect(wdpath, "cuizaixu_lab")){
  interfileFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/interdataFolder_ABCD'
  functionFolder <- '/ibmgpfs/cuizaixu_lab/xuxiaoyu/SC_development/Rcode_SCdevelopment/gamfunction'
}else{
  interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_ABCD'
  functionFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
  FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_ABCD_final/SA12'
  demopath<-"D:/xuxiaoyu/open_dataset_information/ABCD/info"
}

source(paste0(functionFolder, '/SCrankcorr.R'))
source("D:/xuxiaoyu/DMRI_network_development/Normative_model/functions/plotmatrix.R")
# input consistency threshold used to filter spurious streamlines
CVthr=75
SCdata.sum.merge<-readRDS(paste0(interfileFolder, "/SCdata.diw_CV", CVthr, ".rds"))
SCdata <- SCdata.sum.merge %>% drop_na(c(2:(1+elementnum)))
SCdata$eventname <- str_split_i(SCdata$scanID, "_ses-", 2)
# Clean data, keep participants with over one measurements.
SCdata <- SCdata %>% group_by(subID) %>% filter(n() > 1)
SCdata <- SCdata %>% group_by(subID) %>% mutate(
  age_bp = age[eventname=="baselineYear1Arm1"],
  age_wp = age - age_bp
)

summary(SCdata[,c("age_bp", "age_wp")])
# age_bp           age_wp     
# Min.   : 8.917   Min.   :0.000  
# 1st Qu.: 9.333   1st Qu.:0.000  
# Median : 9.917   Median :0.750  
# Mean   : 9.912   Mean   :1.022  
# 3rd Qu.:10.417   3rd Qu.:2.000  
# Max.   :11.000   Max.   :3.167  

SCdata$subID <- as.factor(SCdata$subID)
SCdata$subID <- droplevels(SCdata$subID)

# Fit model and compute individual slopes
lmmfit <- function(dep.var, dataname, bpagevar, wpagevar, covariates, subvar){
  #Fit the lm
  lm.data <- get(dataname)
  lm.data$tmp <- lm.data[[dep.var]]
  lm.data <- lm.data %>% filter((tmp < mean(tmp)+3*sd(tmp)), (tmp > mean(tmp)-3*sd(tmp)))
  
  if (is.na(covariates)){
    modelformula <- as.formula(sprintf("%1$s ~ 1 + %2$s + %3$s + (1 | %4$s)", dep.var, bpagevar, wpagevar, subvar))
  }else{
    modelformula <- as.formula(sprintf("%1$s ~ 1 + %2$s + %3$s +  %4$s + (1 | %5$s)", dep.var, bpagevar, wpagevar, covariates, subvar))
  }
  
  lm.mod <- lmer(modelformula, data = lm.data, control = lmerControl(calc.derivs = FALSE),REML=TRUE)
  SingularInfo <- isSingular(lm.mod)
  # fixed effect
  fixed=summary(lm.mod)$coefficient
  age_bp_beta <- fixed[2, 1]
  age_wp_beta <- fixed[3, 1]
  
  age_bp_T <- fixed[2, 4]
  age_wp_T <- fixed[3, 4]
  age_bp_P <- fixed[2, 5]
  age_wp_P <- fixed[3, 5]
  
  #rand_int=ranef(lm.mod)[[subvar]]
  #slopedf <- data.frame(subID =rownames(rand_int), slope = rand_int[[wpagevar]]+age_wp_beta)
  #names(slopedf) <- c("subID", paste0("slope.", dep.var))
  df <- data.frame(dep.var, age_bp_beta, age_bp_T, age_wp_beta, age_wp_T, age_bp_P, age_wp_P, SingularInfo)
  #result <- list(df, slopedf)
  
  return(df)
}

# Fit models
dataname <- "SCdata"
bpagevar <- "age_bp"
wpagevar <- "age_wp"
covariates <- "sex + mean_fd"
subvar <- "subID"
stats.df <- list()
slope.list <- list()
for (i in 1:78){
  dep.var <- paste0("SC.", i, "_h")
  result.tmp <- lmmfit(dep.var, dataname, bpagevar, wpagevar, covariates, subvar)
  stats.df[[i]] <- result.tmp
}

stats.df <- do.call(rbind, stats.df)
summary(stats.df)

SCrankcorr(stats.df, "age_wp_T", 12)
# ds.resolution Interest.var r.spearman   p.spearman
# 1            12     age_wp_T  -0.5661642 6.580146e-08


stats.df$age_wp_P.fdr <- p.adjust(stats.df$age_wp_P, method = "fdr")
sum(stats.df$age_wp_P.fdr < 0.05) #70
saveRDS(stats.df, paste0(interfileFolder, "/LMM_withinage_effect_SA12.rds"))

fig.agewp <- plotmatrix(dataname="stats.df", variable="age_wp_T", ds.resolution=12, Pvar="age_wp_P.fdr", NAcol="white", lmthr=NA, 
                        axeslabels=NULL, axeslabelsGap=T, linerange_frame=NA, PaletteSet=NA, Pvar.noFDR=NA)
fig.agewp
filename<-paste0(FigureFolder,'/CV', CVthr, "/Matrix12_sumSCinvnode_gamstats_Age8_22/AgeWP_Tscore_delLM_CV75_siteall.tiff")
ggsave(filename, fig.agewp,  height = 16, width = 18, units = "cm")

## 2) Fit trajectories by cognition and p-factor
#################################
# load data
SA12_10 <- read.csv(paste0(interfileFolder, "/SA12_10.csv"))
plotdata <- readRDS(paste0(interfileFolder, '/plotdatasum.df_SA12_sumSCinvnode_siteall_CV', CVthr,'.rds'))
SCdata<-readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatgam_age_sex_meanfd.rds'))
source(paste0(functionFolder, "/gamminteraction_linear.R"))

Behavior <- read.csv(paste0(demopath, '/DemodfScreenFinal.csv'))
Behavior$siteID <- gsub("site0", "site", Behavior$siteID)
summary(Behavior[,c("sex", "eventname", "handness")])
factorvar<-c("sex", "eventname", "handness")
for (i in 1:3){
  Behavior[, factorvar[i]] <-as.factor(Behavior[, factorvar[i]])
}
Behavior <- Behavior %>% filter(eventname=="baseline_year_1_arm_1")

SCdata <- SCdata %>% left_join(select(Behavior, subID, nihtbx_fluidcomp_uncorrected, GENERAL), by="subID")

SCdata.diw <- SCdata
for (x in 1:78){
  region <- grep("SC.", names(SCdata), value = T)[x]
  plotdata.tmp <- plotdata[plotdata$SC_label==paste0("SC.", x, "_h"), ]
  SCstrength.diw <- SCdata[,region] / plotdata.tmp$fit[1]
  SCdata.diw[,region] <- SCstrength.diw
}
SCdata.diw[,grep("SC.", names(SCdata), value = T)] <- lapply(SCdata.diw[,grep("SC.", names(SCdata), value = T)], as.numeric)
SCdata.diw$sex <- as.factor(SCdata.diw$sex)
SCdata.diw$eventname <- str_split_i(SCdata.diw$scanID, "_ses-", 2)
SCdata.diw <- SCdata.diw %>% group_by(subID) %>% filter(n() > 1)
SCdata.diw <- SCdata.diw %>% group_by(subID) %>% mutate(
  age_bp = age[eventname=="baselineYear1Arm1"],
  age_wp = age - age_bp
)

summary(SCdata.diw[,c("age_bp", "age_wp")])
# age_bp           age_wp     
# Min.   : 8.917   Min.   :0.000  
# 1st Qu.: 9.333   1st Qu.:0.000  
# Median : 9.917   Median :0.750  
# Mean   : 9.912   Mean   :1.022  
# 3rd Qu.:10.417   3rd Qu.:2.000  
# Max.   :11.000   Max.   :3.167  

SCdata.diw$subID <- as.factor(SCdata.diw$subID)
SCdata.diw$subID <- droplevels(SCdata.diw$subID)

## Cognition
##########################
dataname <- "SCdata.diw"
smooth_var <- "age_wp"
int_var <- "nihtbx_fluidcomp_uncorrected"
covariates <- "age_bp + sex+mean_fd"
subvar <- "subID"
increments = 1000
stats.all <- list(); plot.all <- list()
for (i in 1:78){
  
  region <- paste0("SC.", i, "_h")
  # low
  int_var.predict.percentile = 0.1
  result.low <- lmm.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                int_var.predict.percentile, 
                                                covariates, subvar, increments, stats_only=F, 
                                                plotdf=NA)
  statsresult <- result.low[[1]]
  stats.all[[i]] <- statsresult
  plot.low <- result.low[[2]]
  plot.low$label <- region
  plot.low$cognitionlevel <- "low"
  # high
  int_var.predict.percentile = 0.9
  result.high <- lmm.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                 int_var.predict.percentile, 
                                                 covariates, subvar, increments, stats_only=F, 
                                                 plotdf=NA)
  plot.high <- result.high[[2]]
  plot.high$label <- region
  plot.high$cognitionlevel <- "high"
  ## Merge
  plotdf <- rbind(plot.low, plot.high)
  plot.all[[i]] <- plotdf
}

saveRDS(plot.all, paste0(interfileFolder, '/plotdata_high90_low10_', int_var, '_develop_WithinPerson.rds'))

# merge plotdf
plotdf <- do.call(rbind, plot.all)
plotdf <- merge(plotdf, SA12_10, by.x = "label", by.y = "SC_label")
plotdf.decile.low <- plotdf %>% filter(cognitionlevel=="low") %>% 
  group_by(decile, age_wp) %>%
  summarise(fit.avg = mean(fit), decile=mean(decile))
plotdf.decile.low$cognitionlevel <- "low"

plotdf.decile.high <- plotdf %>% filter(cognitionlevel=="high") %>% 
  group_by(decile, age_wp) %>%
  summarise(fit.avg = mean(fit), decile=mean(decile))
plotdf.decile.high$cognitionlevel <- "high"
plotdf.decile <- rbind(plotdf.decile.low, plotdf.decile.high)
# > summary(plotdf.decile$fit.avg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8893  0.9703  1.0026  1.0003  1.0262  1.1040 
# Save out plots
colorid <- rev(brewer.pal(10, "RdBu"))
for (i in 1:10){
  plotdf.tmp <- plotdf.decile[plotdf.decile$decile==i,]
  colorindex <- colorid[i]
  if (i == 1){
    mytheme=theme(axis.text=element_text(size=21, color="black"), 
                  axis.title =element_text(size=21),aspect.ratio = 1,
                  axis.line = element_line(linewidth=0.5),axis.ticks = element_line(linewidth=0.5),
                  plot.background=element_rect(fill="transparent"),
                  panel.border = element_rect(fill=NA, color="transparent"),
                  panel.background=element_rect(fill="transparent", color="transparent"), legend.position = "none")
    
  }else{
    mytheme=theme(axis.text.x =element_text(size=21, color="black"), 
                  axis.text.y =element_text(size=21, color="transparent"),
                  axis.title.x =element_text(size=21),
                  axis.title.y =element_text(size=21, colour = "transparent"),
                  aspect.ratio = 1,
                  axis.line.x = element_line(linewidth=0.5), 
                  axis.line.y = element_line(linewidth=0.5, colour = "transparent"),
                  axis.ticks.x = element_line(linewidth=0.5),
                  axis.ticks.y = element_line(linewidth=0.5, colour = "transparent"),
                  plot.background=element_rect(fill="transparent"),
                  panel.grid=element_line(linewidth=0.5, colour = "transparent"),
                  panel.border = element_rect(fill=NA, color="transparent"),
                  panel.background=element_rect(fill="transparent", color="transparent"), legend.position = "none")
    
  }
  
  
  Fig <- ggplot(data=plotdf.tmp)+
    geom_line(aes(x=age_wp, y=fit.avg, group=cognitionlevel, linetype=cognitionlevel), 
              linewidth=1.2, color=colorindex)+
    scale_linetype_manual(values=c("solid", "dashed"))+
    scale_y_continuous(breaks = c(0.9, 1.0, 1.1), limits=c(0.88, 1.11))+
    labs(x=NULL, y="SC strength (ratio)", color="Cognition")+
    mytheme
  
  print(Fig)
  ggsave(paste0(FigureFolder, '/CV75/cognition/', int_var, '_base/Interaction_linear/developmentcurve_decile', i, '.tiff'),Fig, width = 10, height = 10, units = "cm")
  ggsave(paste0(FigureFolder, '/CV75/cognition/', int_var, '_base/Interaction_linear/developmentcurve_decile', i, '.svg'),Fig, width = 8.5, height = 9.5, units = "cm")
  
}

statsresult.df <- do.call(rbind, stats.all)
statsresult.df <- as.data.frame(statsresult.df)
statsresult.df$T.int <- as.numeric(statsresult.df$T.int)
SCrankcorr(statsresult.df, "T.int", 12)
#########################

## P-factor
##########################
dataname <- "SCdata.diw"
smooth_var <- "age_wp"
int_var <- "GENERAL"
covariates <- "age_bp + sex+mean_fd"
subvar <- "subID"
increments = 1000
stats.all <- list(); plot.all <- list()
for (i in 1:78){
  
  region <- paste0("SC.", i, "_h")
  # low
  int_var.predict.percentile = 0.1
  result.low <- lmm.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                 int_var.predict.percentile, 
                                                 covariates, subvar, increments, stats_only=F, 
                                                 plotdf=NA)
  statsresult <- result.low[[1]]
  stats.all[[i]] <- statsresult
  plot.low <- result.low[[2]]
  plot.low$label <- region
  plot.low$pfactorlevel <- "low"
  # high
  int_var.predict.percentile = 0.9
  result.high <- lmm.predict.covariateinteraction(region, dataname, smooth_var, int_var,
                                                  int_var.predict.percentile, 
                                                  covariates, subvar, increments, stats_only=F, 
                                                  plotdf=NA)
  plot.high <- result.high[[2]]
  plot.high$label <- region
  plot.high$pfactorlevel <- "high"
  ## Merge
  plotdf <- rbind(plot.low, plot.high)
  plot.all[[i]] <- plotdf
}

saveRDS(plot.all, paste0(interfileFolder, '/plotdata_high90_low10_', int_var, '_develop_WithinPerson.rds'))

# merge plotdf
plotdf <- do.call(rbind, plot.all)
plotdf <- merge(plotdf, SA12_10, by.x = "label", by.y = "SC_label")
plotdf.decile.low <- plotdf %>% filter(pfactorlevel=="low") %>% 
  group_by(decile, age_wp) %>%
  summarise(fit.avg = mean(fit), decile=mean(decile))
plotdf.decile.low$pfactorlevel <- "low"

plotdf.decile.high <- plotdf %>% filter(pfactorlevel=="high") %>% 
  group_by(decile, age_wp) %>%
  summarise(fit.avg = mean(fit), decile=mean(decile))
plotdf.decile.high$pfactorlevel <- "high"
plotdf.decile <- rbind(plotdf.decile.low, plotdf.decile.high)
# > summary(plotdf.decile$fit.avg)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9032  0.9769  1.0038  0.9992  1.0216  1.0802 
# Save out plots
colorid <- rev(brewer.pal(10, "RdBu"))
for (i in 1:10){
  plotdf.tmp <- plotdf.decile[plotdf.decile$decile==i,]
  colorindex <- colorid[i]
  if (i == 1){
    mytheme=theme(axis.text=element_text(size=21, color="black"), 
                  axis.title =element_text(size=21),aspect.ratio = 1,
                  axis.line = element_line(linewidth=0.5),axis.ticks = element_line(linewidth=0.5),
                  plot.background=element_rect(fill="transparent"),
                  panel.border = element_rect(fill=NA, color="transparent"),
                  panel.background=element_rect(fill="transparent", color="transparent"), legend.position = "none")
    
  }else{
    mytheme=theme(axis.text.x =element_text(size=21, color="black"), 
                  axis.text.y =element_text(size=21, color="transparent"),
                  axis.title.x =element_text(size=21),
                  axis.title.y =element_text(size=21, colour = "transparent"),
                  aspect.ratio = 1,
                  axis.line.x = element_line(linewidth=0.5), 
                  axis.line.y = element_line(linewidth=0.5, colour = "transparent"),
                  axis.ticks.x = element_line(linewidth=0.5),
                  axis.ticks.y = element_line(linewidth=0.5, colour = "transparent"),
                  plot.background=element_rect(fill="transparent"),
                  panel.grid=element_line(linewidth=0.5, colour = "transparent"),
                  panel.border = element_rect(fill=NA, color="transparent"),
                  panel.background=element_rect(fill="transparent", color="transparent"), legend.position = "none")
    
  }
  
  
  Fig <- ggplot(data=plotdf.tmp)+
    geom_line(aes(x=age_wp, y=fit.avg, group=pfactorlevel, linetype=pfactorlevel), 
              linewidth=1.2, color=colorindex)+
    scale_linetype_manual(values=c("dashed", "solid"))+
    scale_y_continuous(breaks = c(0.9, 1.0, 1.1), limits=c(0.88, 1.11))+
    labs(x=NULL, y="SC strength (ratio)", color="p-factor")+
    mytheme
  
  print(Fig)
  ggsave(paste0(FigureFolder, '/CV75/Disease/pFactor/Interaction_linear/developmentcurve_decile', i, '.tiff'),
         Fig, width = 10, height = 10, units = "cm")

}

statsresult.df <- do.call(rbind, stats.all)
statsresult.df <- as.data.frame(statsresult.df)


