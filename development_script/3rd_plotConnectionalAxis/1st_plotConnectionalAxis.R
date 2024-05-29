library(ggplot2)
library(pals)
library(scales)
library(grDevices)
FigureFolder<-'/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Figure_HCPD_final/SA12'

# Connectional rank
Matrix.tmp <- matrix(NA, nrow = 12, ncol=12)
for (x in 1:12){
  for (y in 1:12){
    Matrix.tmp[x,y] <- (x+y)^2+(x-y)^2
  }
}
SCrank <- rank(Matrix.tmp[lower.tri(Matrix.tmp, diag = T)], ties.method = "average")
Matrix.tmp[lower.tri(Matrix.tmp, diag = T)] <- SCrank
Matrix.tmp[upper.tri(Matrix.tmp)] <- t(Matrix.tmp)[upper.tri(Matrix.tmp)]
colnames(Matrix.tmp) <-seq(1, 12)
rownames(Matrix.tmp) <-seq(1, 12)
matrixtmp.df <- as.data.frame(Matrix.tmp)
matrixtmp.df$nodeid <- seq(1, 12)
matrixtmp.df.melt <- melt(matrixtmp.df,id.vars=c("nodeid"))
matrixtmp.df.melt$variable<-as.numeric(matrixtmp.df.melt$variable)
matrixtmp.df.melt$nodeid<-0-matrixtmp.df.melt$nodeid
matrixtmp.df.melt$value<-as.numeric(matrixtmp.df.melt$value)
linerange_frame<-data.frame(x=c(0.5,12+0.5), ymin =rep(-12-0.5, times=2), ymax =rep(-0.5, times=2),
                            y=c(-0.5, -12-0.5), xmin=rep(0.5, times=2), xmax=rep(12+0.5, times=2))

ggplot(data =matrixtmp.df.melt)+
  geom_tile(aes(x=variable, y=nodeid, fill = value, color=value))+
  scale_fill_distiller(type="seq", palette = "RdBu",na.value = "grey")+
  scale_color_distiller(type="seq", palette = "RdBu",na.value = "grey")+
  geom_linerange(data=linerange_frame, aes(y=y, xmin =xmin, xmax =xmax), color="black", linewidth=0.5)+
  geom_linerange(data=linerange_frame, aes(x=x, ymin =ymin, ymax =ymax), color="black", linewidth=0.5)+
  geom_segment(aes(x = 0.5 , y = -0.5 , xend = 12+0.5 ,yend = -12-0.5), color="black", linewidth=0.5)+
  ggtitle(label = "Connectional rank")+labs(x=NULL, y=NULL)+
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
filename<-paste0(FigureFolder, "/Connectionrank_squareDistance.tiff")
ggsave(filename, height = 18, width = 20, units = "cm")

## save out RGB values of colormap pals::ocean.matter
colors <- ocean.matter(376)
rgb_values <- t(col2rgb(colors))
rgb_values_normalized <- rgb_values / 255
rgb_df <- as.data.frame(rgb_values_normalized)
colnames(rgb_df) <- c("R", "G", "B")
rgb_df$U <- 1
write.csv(rgb_df, paste0(interfileFolder, "/oceanmatter376_rgb.csv"), row.names = FALSE)

## purple --> yellow
low_color <- "goldenrod1"
mid_color <- "white"
high_color <- "#6f1282"

palette_func <- colorRampPalette(c(low_color, mid_color, high_color))
colors <- palette_func(376)
rgb_values <- t(col2rgb(colors))
rgb_values_normalized <- rgb_values / 255
rgb_df <- as.data.frame(rgb_values_normalized)
colnames(rgb_df) <- c("R", "G", "B")
rgb_df$U <- 1

write.csv(rgb_df, paste0(interfileFolder, "/goldenrod_white376_rgb.csv"), row.names = FALSE)

color_df <- data.frame(
  x = 1:376,
  color = colors
)
ggplot(color_df, aes(x = x, y = 1, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Custom Color Palette",
       x = NULL,
       y = NULL)


