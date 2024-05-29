# compute Euclidean distance between nodes
rm(list=ls())
SA12.coords.lh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA12_lh.csv", header=F) #coordinates of SA12 parcel centroids on the freesurfer sphere
SA12.coords.rh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA12_rh.csv", header=F) #coordinates of SA12 parcel centroids on the freesurfer sphere
distance.lh <- matrix(NA, 12, 12)
for (x in 1:12){
  for (y in 1:12){
    distance.tmp = ((SA12.coords.lh$V1[x]-SA12.coords.lh$V1[y])^2+
      (SA12.coords.lh$V2[x]-SA12.coords.lh$V2[y])^2+
      (SA12.coords.lh$V3[x]-SA12.coords.lh$V3[y])^2)^0.5
    distance.lh[x,y] = distance.tmp+1
  }
}
distance.rh <- matrix(NA, 12, 12)
for (x in 1:12){
  for (y in 1:12){
    distance.tmp = ((SA12.coords.rh$V1[x]-SA12.coords.rh$V1[y])^2+
                      (SA12.coords.rh$V2[x]-SA12.coords.rh$V2[y])^2+
                      (SA12.coords.rh$V3[x]-SA12.coords.rh$V3[y])^2)^0.5
    distance.rh[x,y] = distance.tmp+1
  }
}

image(distance.lh)
image(distance.rh)
distance.avg <- (distance.lh+distance.rh)/2
image(distance.avg)
distance.avg.df <- data.frame(Edistance=distance.avg[lower.tri(distance.avg, diag=T)])
distance.avg.df$SC_label <- paste0("SC.", 1:78, "_h")
write.csv(distance.avg.df, "/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_HCPD/average_EuclideanDistance_12.csv", row.names = F)
write.csv(distance.avg.df, "/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/interdataFolder_ABCD/average_EuclideanDistance_12.csv", row.names = F)
