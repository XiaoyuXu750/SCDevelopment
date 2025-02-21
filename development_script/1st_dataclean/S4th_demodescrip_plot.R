# generate age distribution plots and demographic tables.
library(ggplot2)
library(RColorBrewer)
library('reshape2')
library(gdata)
library(tableone)
display.brewer.all()

rm(list = ls())
interfileFolderABCD <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD'
interfileFolderHCPD <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD'
FigureFolderABCD <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_ABCD_final'
FigureFolderHCPD <- '/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final'
dataHCPD <- readRDS(paste0(interfileFolderHCPD, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))
dataABCD <- readRDS(paste0(interfileFolderABCD, "/SCdata_SA12_CV75_sumSCinvnode.sum.msmtcsd.combatage.rds"))
dataABCD$age <- dataABCD$age / 12
# HCPD
ggplot(data = dataHCPD, aes(age, y = ..count..)) +
  geom_histogram(binwidth = 1, color = "black", fill = "#B4D3E7", boundary = 5, position = position_dodge(width = 0.8)) +
  labs(x = "Age (years)", y = "Frequency", title = paste0("HCPD, N=", nrow(dataHCPD))) +
  scale_x_continuous(limits = c(8, 23), breaks = c(8,11,14,17,20,23)) +
  #scale_y_continuous(limits = c(0, 70), breaks = c(0,10, 20,30, 40,50,60,70)) +
  #geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  labs(y=NULL)+theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
    plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.7,
    plot.title = element_text(color = "black", size = 14, hjust = 0.5),
    axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(color = "black", size = 14),
    legend.position = "none")

ggsave(paste(FigureFolderHCPD, '/Fig1_Age_distribution_count_HCPD.tiff', sep = ''), width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderHCPD, '/Fig1_Age_distribution_count_HCPD.svg', sep = ''), width = 12, height = 10, units = "cm")

# ABCD
fillcolor = c("#83B7D7", "#B4D3E7")
ggplot(data = dataABCD, aes(age, y = ..count.., fill = eventname)) +
  geom_histogram(binwidth = 0.5, color = "black", boundary = 5, position = "stack",linewidth=0.5) +
  labs(x = "Age (years)", y = NULL, title = "ABCD, N=7,104") +
  #scale_x_continuous(limits = c(9, 16), breaks = c(9,10,11,12,13,14,15,16)) +
  #scale_y_continuous(limits = c(0, 50), breaks = c(0,10, 20,30,40, 50)) +
  scale_fill_manual(values = fillcolor)+
  #geom_hline(aes(yintercept = 10), colour = "red", linetype="dashed")+
  theme_classic()+
  theme(panel.background = element_rect(fill="transparent"),
        plot.background = element_rect(fill = "transparent",colour = NA),aspect.ratio = 0.75,
        plot.title = element_text(color = "black", size = 14, hjust = 0.5),
        axis.title = element_text(color = "black", size = 14),axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(color = "black", size = 14),
        legend.position = "none")

ggsave(paste(FigureFolderABCD, '/Fig1_Age_distribution_count_ABCD.tiff', sep = ''),  width = 12, height = 12, units = "cm")
ggsave(paste(FigureFolderABCD, '/Fig1_Age_distribution_count_ABCD.svg', sep = ''), width = 12, height = 10, units = "cm")

length(which(table(dataABCD$subID)==2)) # 2570 participants have 2 visits

## Description for demographic information
dataHCPD <- within(dataHCPD, {
  nih_fluidcogcomp_unadjusted[nih_fluidcogcomp_unadjusted<mean(nih_fluidcogcomp_unadjusted, na.rm=T)-3*sd(nih_fluidcogcomp_unadjusted, na.rm=T) | nih_fluidcogcomp_unadjusted>mean(nih_fluidcogcomp_unadjusted, na.rm=T)+3*sd(nih_fluidcogcomp_unadjusted, na.rm=T)] <- NA
})
summary(dataHCPD$nih_fluidcogcomp_unadjusted)
HCPD_vars <- c("age", "sex", "handnessfactor", "race_ethnicity", "mean_fd", "nih_fluidcogcomp_unadjusted", "site", "ICV", "income.adj")
HCPD_tableone <- CreateTableOne(HCPD_vars, data = dataHCPD, factorVars=c("sex", "handnessfactor", "race_ethnicity", "site"))
HCPD_tableone <- as.data.frame(print(HCPD_tableone))
write.csv(HCPD_tableone, paste0(FigureFolderHCPD, '/demographic_info.csv'))

ABCD_vars <- c("age", "sex","eventname", "handness", "race_ethnicity", "mean_fd", "nihtbx_fluidcomp_uncorrected","GENERAL", "siteID", "smri_vol_scs_intracranialv", "income.adj")
ABCD_tableone <- CreateTableOne(ABCD_vars,strata="eventname", data = dataABCD, factorVars=c("sex", "handness", "race_ethnicity", "siteID"), test=F)
ABCD_tableone <- as.data.frame(print(ABCD_tableone))
write.csv(ABCD_tableone, paste0(FigureFolderABCD, '/demographic_info.csv'))



