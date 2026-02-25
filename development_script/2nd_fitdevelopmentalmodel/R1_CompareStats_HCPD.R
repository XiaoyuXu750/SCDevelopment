library(effectsize)
library(psych)
library(tidyverse)
library(parallel)
library(scales)
library(openxlsx)
library(gratia)
library(ggplot2)
library(RColorBrewer)
rm(list = ls())
# Set path
resultFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/results_HCPD'
interfileFolder <- 'D:/xuxiaoyu/DMRI_network_development/SC_development/interdataFolder_HCPD'
functionFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction'
FigureFolder<-'D:/xuxiaoyu/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'
# Load data
CVthr = 75
gamresult <- readRDS(paste0(interfileFolder, '/gamresults78_sumSCinvnode_over8_CV', CVthr,'_scale_TRUE.rds'))
SCdata <- readRDS(paste0(interfileFolder, '/SCdata_SA12_CV', CVthr,'_sumSCinvnode.sum.msmtcsd.combatgam.rds'))

# 1. Direct compare and compute cohen's d
matindex <- matrix(NA, 12, 12)
matindex[lower.tri(matindex, diag = T)] <- 1:78
Type.S <- as.numeric(matindex[1:4, 1:4])
Type.SA <- as.numeric(matindex[5:12, 1:4])
Type.A <- as.numeric(matindex[5:12, 5:12])

Type.S <- Type.S[! is.na(Type.S)]
Type.SA <- Type.SA[! is.na(Type.SA)]
Type.A <- Type.A[! is.na(Type.A)]

gamresult$Type <- NA
gamresult$Type[Type.S] <- "SS"
gamresult$Type[Type.A] <- "AA"
gamresult$Type[Type.SA] <- "SA"

# ANOVA test for three type
## Partial Rsq
##################################################
anova_mod <- aov(partialRsq ~ Type, data = gamresult)
summary(anova_mod)
#            Df  Sum Sq  Mean Sq F value   Pr(>F)    
# Type         2 0.01758 0.008790   15.47 2.38e-06 ***
# Residuals   75 0.04263 0.000568          
# post-hoc
posthoc <- TukeyHSD(anova_mod)

resultdf <- data.frame(Pairs=c("SS_AA", "SS_SA", "SA_AA"), Tvalue=rep(NA, 3),df=rep(NA, 3), Pvalue=rep(NA, 3), cohend=rep(NA, 3))
# SS VS. AA
df.tmp <- gamresult %>% filter(Type %in% c("SS", "AA"))
ttest_SS_AA <- t.test(partialRsq ~ Type, data=df.tmp, paired = F)
resultdf$Tvalue[1] <- ttest_SS_AA$statistic; resultdf$Pvalue[1] <- ttest_SS_AA$p.value; resultdf$df[1] <- ttest_SS_AA$parameter
cohend.tmp <- cohens_d(partialRsq ~ Type, data=df.tmp)
resultdf$cohend[1] <- cohend.tmp$Cohens_d
# SS VS. SA
df.tmp <- gamresult %>% filter(Type %in% c("SS", "SA"))
ttest_SS_SA <- t.test(partialRsq ~ Type, data=df.tmp, paired = F)
resultdf$Tvalue[2] <- ttest_SS_SA$statistic; resultdf$Pvalue[2] <- ttest_SS_SA$p.value; resultdf$df[2] <- ttest_SS_SA$parameter
cohend.tmp <- cohens_d(partialRsq ~ Type, data=df.tmp)
resultdf$cohend[2] <- cohend.tmp$Cohens_d
# SA VS. AA
df.tmp <- gamresult %>% filter(Type %in% c("SA", "AA"))
ttest_SA_AA <- t.test(partialRsq ~ Type, data=df.tmp, paired = F)
resultdf$Tvalue[3] <- ttest_SA_AA$statistic; resultdf$Pvalue[3] <- ttest_SA_AA$p.value; resultdf$df[3] <- ttest_SA_AA$parameter
cohend.tmp <- cohens_d(partialRsq ~ Type, data=df.tmp)
resultdf$cohend[3] <- cohend.tmp$Cohens_d

resultdf$Pvalue.adj <- p.adjust(resultdf$Pvalue, method="fdr")

#   Pairs    Tvalue       df       Pvalue     cohend   Pvalue.adj
# 1 SS_AA -3.108609 12.65763 0.0085438207 -1.2366230 0.0085438207
# 2 SS_SA -4.855129 13.40250 0.0002882196 -1.9100801 0.0008646588
# 3 SA_AA  3.067399 64.52086 0.0031542442  0.7466638 0.0047313663

# Plot bar
stats_df <- gamresult %>% group_by(Type) %>%
  summarise(
    mean_Rsq = mean(partialRsq, na.rm = TRUE),
    sd_Rsq = sd(partialRsq, na.rm = TRUE),
    .groups = "drop"
  )

stats_df$Type <- factor(stats_df$Type, levels = c("SS", "SA", "AA"))
barplot <- ggplot(stats_df, aes(x = Type, y = mean_Rsq, fill = Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_Rsq - sd_Rsq, ymax = mean_Rsq + sd_Rsq),
    position = position_dodge(width = 0.8),
    width = 0.25,
    color = "black"
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)))+
  labs(
    x = "Type of connections",
    y = expression("Age effects (partial "*italic(R)^2 * ")"),
    fill = NULL
  ) +
  scale_fill_manual(values = c("#2166AC", "#F4A582", "#B2182B")) +
  theme_classic()+
  theme(axis.text=element_text(size=18.4, color="black"), 
        axis.title =element_text(size=18.4, color="black"),aspect.ratio = 1,
        plot.background=element_rect(fill="white"),
        panel.background=element_rect(fill="white"),
        plot.title = element_text(size=16, hjust = 0.5),
        legend.text = element_text(size=18.5, color="black"), legend.title = element_text(size=20, color="black"))

barplot
ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/PartialRsq_barplot.tiff"), barplot, width=17, height =15.7, units = "cm")


# 2. Permutation test
pairlist <- list(c("SS", "AA"), c("SS", "SA"), c("SA", "AA"))
permutationresult <- list()
for (p in 1:1000){
  set.seed(925 + p)
  gamresult.perm <- gamresult
  gamresult.perm$partialRsq <- gamresult.perm$partialRsq[sample(1:78, 78)]
  
  resultdf.null <- data.frame(Pairs=c("SS_AA", "SS_SA", "SA_AA"), Tvalue=rep(NA, 3),df=rep(NA, 3), Pvalue=rep(NA, 3), cohend=rep(NA, 3))
  for (t in 1:3){
    df.tmp <- gamresult.perm %>% filter(Type %in% pairlist[[t]])
    ttest_tmp <- t.test(partialRsq ~ Type, data=df.tmp, paired = F)
    resultdf.null$Tvalue[t] <- ttest_tmp$statistic; resultdf.null$Pvalue[t] <- ttest_tmp$p.value; resultdf.null$df[t] <- ttest_tmp$parameter
    cohend.tmp <- cohens_d(partialRsq ~ Type, data=df.tmp)
    resultdf.null$cohend[t] <- cohend.tmp$Cohens_d
  }
 
  permutationresult[[p]] <- resultdf.null
}

# Null distribution
permutationresult.df <- do.call(rbind, lapply(permutationresult, function(x){
  df <- data.frame(SS_AA=x$cohend[1], SS_SA=x$cohend[2], SA_AA=x$cohend[3])
  return(df)
}))

# Pvalues
p.SS_AA <- sum(resultdf$cohend[1] >= permutationresult.df$SS_AA) / 1000 # 0.002
p.SS_SA <- sum(resultdf$cohend[2] >= permutationresult.df$SS_SA) / 1000 # 0.00
p.SA_AA <- sum(resultdf$cohend[3] <= permutationresult.df$SA_AA) / 1000 # 0.002

# SS VS AA
ggplot(data = permutationresult.df, aes(SS_AA, y = ..count..)) +
  geom_histogram(color = "black", fill = "#B4D3E7", position = position_dodge(width = 0.8)) +
  labs(x = "Null distribution", y = "Frequency", title = "SS VS. AA") +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  geom_vline(aes(xintercept = resultdf$cohend[1]), colour = "red", linetype="solid", linewidth=1)+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1,
        plot.title = element_text(color = "black", size = 13.5, hjust = 0.5),
        axis.title = element_text(color = "black", size = 13.5),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 13.5),
        legend.position = "none")

ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/SS_AA_permutation.tiff"), width = 12, height = 10, units = "cm")

# SS VS SA
ggplot(data = permutationresult.df, aes(SS_SA, y = ..count..)) +
  geom_histogram(color = "black", fill = "#B4D3E7", position = position_dodge(width = 0.8)) +
  labs(x = "Null distribution", y = "Frequency", title = "SS VS. SA") +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  geom_vline(aes(xintercept = resultdf$cohend[2]), colour = "red", linetype="solid", linewidth=1)+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1,
        plot.title = element_text(color = "black", size = 13.5, hjust = 0.5),
        axis.title = element_text(color = "black", size = 13.5),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 13.5),
        legend.position = "none")

ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/SS_SA_permutation.tiff"), width = 12, height = 10, units = "cm")

# SA VS AA
ggplot(data = permutationresult.df, aes(SA_AA, y = ..count..)) +
  geom_histogram(color = "black", fill = "#B4D3E7", position = position_dodge(width = 0.8)) +
  labs(x = "Null distribution", y = "Frequency", title = "SA VS. AA") +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  geom_vline(aes(xintercept = resultdf$cohend[3]), colour = "red", linetype="solid", linewidth=1)+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1,
        plot.title = element_text(color = "black", size = 13.5, hjust = 0.5),
        axis.title = element_text(color = "black", size = 13.5),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 13.5),
        legend.position = "none")

ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/SA_AA_permutation.tiff"), width = 12, height = 10, units = "cm")
######################################################################################

## Second Derivative
#############################################################

## 2.1 Compare with Zero
resultdf.zero <- data.frame(Pairs=c("SS", "SA", "AA"), Tvalue=rep(NA, 3),df=rep(NA, 3), Pvalue=rep(NA, 3), cohend=rep(NA, 3))
# SS VS. 0
df.tmp <- gamresult %>% filter(Type %in% c("SS"))
ttest_SS <- t.test(df.tmp$meanderv2, mu=0)
resultdf.zero$Tvalue[1] <- ttest_SS$statistic; resultdf.zero$Pvalue[1] <- ttest_SS$p.value; resultdf.zero$df[1] <- ttest_SS$parameter
cohend.tmp <- cohens_d(df.tmp$meanderv2, mu=0)
resultdf.zero$cohend[1] <- cohend.tmp$Cohens_d
# SA VS. 0
df.tmp <- gamresult %>% filter(Type %in% c("SA"))
ttest_SA <- t.test(df.tmp$meanderv2, mu=0)
resultdf.zero$Tvalue[2] <- ttest_SA$statistic; resultdf.zero$Pvalue[2] <- ttest_SA$p.value; resultdf.zero$df[2] <- ttest_SA$parameter
cohend.tmp <- cohens_d(df.tmp$meanderv2, mu=0)
resultdf.zero$cohend[2] <- cohend.tmp$Cohens_d
# AA VS. 0
df.tmp <- gamresult %>% filter(Type %in% c("AA"))
ttest_AA <- t.test(df.tmp$meanderv2, mu=0)
resultdf.zero$Tvalue[3] <- ttest_AA$statistic; resultdf.zero$Pvalue[3] <- ttest_AA$p.value; resultdf.zero$df[3] <- ttest_AA$parameter
cohend.tmp <- cohens_d(df.tmp$meanderv2, mu=0)
resultdf.zero$cohend[3] <- cohend.tmp$Cohens_d

resultdf.zero$Pvalue.adj <- p.adjust(resultdf.zero$Pvalue, method="fdr")


# Plot bar
stats_df <- gamresult %>% group_by(Type) %>%
  summarise(
    mean_2ndderv = mean(meanderv2, na.rm = TRUE),
    sd_2ndderv = sd(meanderv2, na.rm = TRUE),
    .groups = "drop"
  )

stats_df$Type <- factor(stats_df$Type, levels = c("SS", "SA", "AA"))
barplot <- ggplot(stats_df, aes(x = Type, y = mean_2ndderv, fill = Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_2ndderv - sd_2ndderv, ymax = mean_2ndderv + sd_2ndderv),
    position = position_dodge(width = 0.8),
    width = 0.25,
    color = "black"
  ) +
  geom_hline(yintercept = 0, linewidth=0.5, color="black")+
  labs(
    x = "Type of connections",
    y = expression("Second derivative"),
    fill = NULL
  ) +
  scale_fill_manual(values = c("#2166AC", "#F4A582", "#B2182B")) +
  theme_classic()+
  theme(axis.text=element_text(size=18.4, color="black"), 
        axis.title =element_text(size=18.4, color="black"),aspect.ratio = 1,
        axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        plot.background=element_rect(fill="white"),
        panel.background=element_rect(fill="white"),
        plot.title = element_text(size=16, hjust = 0.5),
        legend.text = element_text(size=18.5, color="black"), legend.title = element_text(size=20, color="black"))

ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/SecondDeriv_barplot.tiff"), barplot, width=17, height =15.7, units = "cm")

# 2.2 Compare VS Zero
permutationresult.zero <- list()
for (p in 1:1000){
  set.seed(925 + p)
  gamresult.perm <- gamresult
  gamresult.perm$meanderv2 <- gamresult.perm$meanderv2[sample(1:78, 78)]
  
  resultdf.zero.null <- data.frame(Pairs=c("SS", "SA", "AA"), Tvalue=rep(NA, 3),df=rep(NA, 3), Pvalue=rep(NA, 3), cohend=rep(NA, 3))
  
  for (t in 1:3){
    df.tmp <- gamresult.perm %>% filter(Type %in% resultdf.zero$Pairs[[t]])
    ttest_tmp <- t.test(df.tmp$meanderv2, mu=0)
    resultdf.zero.null$Tvalue[t] <- ttest_tmp$statistic; resultdf.zero.null$Pvalue[t] <- ttest_tmp$p.value; 
    resultdf.zero.null$df[t] <- ttest_tmp$parameter
    cohend.tmp <- cohens_d(df.tmp$meanderv2, mu=0)
    resultdf.zero.null$cohend[t] <- cohend.tmp$Cohens_d
  }
  
  permutationresult.zero[[p]] <- resultdf.zero.null
}

# Null distribution
permutationresult.df <- do.call(rbind, lapply(permutationresult, function(x){
  df <- data.frame(SS_AA=x$cohend[1], SS_SA=x$cohend[2], SA_AA=x$cohend[3])
  return(df)
}))

permutationresult.zero.df <- do.call(rbind, lapply(permutationresult.zero, function(x){
  df <- data.frame(SS=x$cohend[1], SA=x$cohend[2], AA=x$cohend[3])
  return(df)
}))

# Pvalues
p.SS_0 <- sum(resultdf.zero$cohend[1] >= permutationresult.zero.df$SS) / 1000 # 0.00
p.SA_0 <- sum(resultdf.zero$cohend[2] <= permutationresult.zero.df$SA) / 1000 # 0.086
p.AA_0 <- sum(resultdf.zero$cohend[3] <= permutationresult.zero.df$AA) / 1000 # 0.00


### Second derivatives VS Zero
# SS VS Zero
ggplot(data = permutationresult.zero.df, aes(SS, y = ..count..)) +
  geom_histogram(color = "black", fill = "#B4D3E7", position = position_dodge(width = 0.8)) +
  labs(x = "Null distribution", y = "Frequency", title = "SS VS. Zero") +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  geom_vline(aes(xintercept = resultdf.zero$cohend[1]), colour = "red", linetype="solid", linewidth=1)+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1,
        plot.title = element_text(color = "black", size = 13.5, hjust = 0.5),
        axis.title = element_text(color = "black", size = 13.5),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 13.5),
        legend.position = "none")

ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/SS_Zero_permutation_2ndderv.tiff"), width = 12, height = 10, units = "cm")

# SA VS Zero
ggplot(data = permutationresult.zero.df, aes(SA, y = ..count..)) +
  geom_histogram(color = "black", fill = "#B4D3E7", position = position_dodge(width = 0.8)) +
  labs(x = "Null distribution", y = "Frequency", title = "SA VS. Zero") +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  geom_vline(aes(xintercept = resultdf.zero$cohend[2]), colour = "red", linetype="solid", linewidth=1)+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1,
        plot.title = element_text(color = "black", size = 13.5, hjust = 0.5),
        axis.title = element_text(color = "black", size = 13.5),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 13.5),
        legend.position = "none")

ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/SA_Zero_permutation_2ndderv.tiff"), width = 12, height = 10, units = "cm")

# AA VS Zero
ggplot(data = permutationresult.zero.df, aes(AA, y = ..count..)) +
  geom_histogram(color = "black", fill = "#B4D3E7", position = position_dodge(width = 0.8)) +
  labs(x = "Null distribution", y = "Frequency", title = "SA VS. AA") +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  geom_vline(aes(xintercept = resultdf.zero$cohend[3]), colour = "red", linetype="solid", linewidth=1)+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 1,
        plot.title = element_text(color = "black", size = 13.5, hjust = 0.5),
        axis.title = element_text(color = "black", size = 13.5),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 13.5),
        legend.position = "none")

ggsave(paste0(FigureFolder, "/CV", CVthr, "/CompareStats/AA_Zero_permutation_2ndderv.tiff"), width = 12, height = 10, units = "cm")

##########################################################


