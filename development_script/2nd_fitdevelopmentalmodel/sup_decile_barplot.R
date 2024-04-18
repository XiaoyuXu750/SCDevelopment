derivative <- merge(derivative, SA12_10, by.x="label_ID", by.y = "SC_label")
derivative <- within(derivative, {
  decile2 <- NA
  decile2[decile <= 5] = 1
  decile2[decile > 5] = 2
})

derivative.df.decile2 <- derivative %>%
  group_by(decile2, age) %>%
  summarise(sig.avg = mean(significant), SCranktype_order=mean(decile))
tmp <- derivative.df.decile2[derivative.df.decile2$decile2==1,]
tmp$decile2 <- 1.5
tmp$sig.avg <- 0
derivative.df.decile2 <- rbind(derivative.df.decile2, tmp)
derivativeplot <- ggplot(data=derivative.df.decile2)+
  geom_bar(aes(x=age, y=c(rep(1,1000),rep(1,1000), rep(0.5,1000)), fill = sig.avg, group=decile2,
               color=sig.avg),stat = "identity", position = "stack")+
  scale_fill_gradient2(high = "#B2182B", low = "white",midpoint=0, na.value = "white", labels=NULL) +
  scale_color_gradient2(high = "#B2182B", low = "white",midpoint=0, na.value = "white",  labels=NULL) +
  scale_y_continuous(breaks = NULL)+
  ylab(NULL)+xlab(NULL)+
  scale_x_continuous(breaks = NULL)+
  theme_classic()+
  theme(axis.text=element_text(size=20, color='black'),
        axis.title = element_text(size = 20),
        axis.line.y=element_line(size=0),
        axis.line.x=element_line(size=0),
        legend.position = 'none')
derivativeplot

lineplot<-ggplot(data=plotdatasum.df.decile, aes(x=age, y=fit.Z, group=decile, color=decile))+
  geom_line(size=2, alpha=0.8)+
  #paletteer::scale_color_paletteer_c("pals::ocean.matter", direction = -1, values=colorbarvalues.meanderiv2, oob = squish) +
  scale_color_distiller(type="seq", palette = "RdBu", direction = -1)+
  labs(x="Age (years)", y="SC strength (z-score)")+
  #scale_color_manual(values = rev(brewer.pal(6, "RdBu")))+
  theme_classic()+
  theme(axis.text=element_text(size=20, color="black"), 
        axis.title =element_text(size=20, color="black"),aspect.ratio = 0.8,
        plot.background=element_rect(fill="transparent"),legend.position = "none",
        panel.background=element_rect(fill="transparent"),
        plot.title = element_text(size=15, hjust = 0.5))
lineplot

# merge plots
allplots <- list(lineplot,derivativeplot)
mergeplot <- cowplot::plot_grid(rel_heights = c(16,4), plotlist = allplots, 
                                align = "v", axis = "lr",greedy=F, ncol = 1, nrow = 2,
                                hjust = -0.5)
mergeplot


ggsave(paste0(FigureFolder,'/CV',CVthr,  '/SA12_decile_sumSCinvnode_fit/devcurve_SCrank_fit.Z_SCtype10_2.tiff'), width=20, height =14, units = "cm")
ggsave(paste0(FigureFolder,'/CV',CVthr,  '/SA12_decile_sumSCinvnode_fit/devcurve_SCrank_fit.Z_SCtype10_2.svg'), dpi=600, width=18, height =15, units = "cm")




