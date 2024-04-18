##Spin Test Parcel Rotation Matrix
source("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/rotate.parcellation.R")
source("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/perm.sphere.p.R")
SA12.coords.lh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA12_lh.csv", header=F) #coordinates of SA12 parcel centroids on the freesurfer sphere
SA12.coords.rh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA12_rh.csv", header=F) #coordinates of SA12 parcel centroids on the freesurfer sphere
perm.id.full <- rotate.parcellation(coord.l = as.matrix(SA12.coords.lh), coord.r = as.matrix(SA12.coords.rh), nrot = 10000) #rotate the SA12 parcellation 10,000 times on the freesurfer sphere to generate spatial nulls for spin-based permutation significance testing 
saveRDS(perm.id.full, "/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA12_sphericalrotations_N10000.rds")

SA17.coords.lh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA17_lh.csv", header=F) #coordinates of SA17 parcel centroids on the freesurfer sphere
SA17.coords.rh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA17_rh.csv", header=F) #coordinates of SA17 parcel centroids on the freesurfer sphere
perm.id.full <- rotate.parcellation(coord.l = as.matrix(SA17.coords.lh), coord.r = as.matrix(SA17.coords.rh), nrot = 10000) #rotate the SA17 parcellation 10,000 times on the freesurfer sphere to generate spatial nulls for spin-based permutation significance testing 
saveRDS(perm.id.full, "/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA17_sphericalrotations_N10000.rds")

SA7.coords.lh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA7_lh.csv", header=F) #coordinates of SA7 parcel centroids on the freesurfer sphere
SA7.coords.rh <- read.csv("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/sphere_SA7_rh.csv", header=F) #coordinates of SA7 parcel centroids on the freesurfer sphere
perm.id.full <- rotate.parcellation(coord.l = as.matrix(SA7.coords.lh), coord.r = as.matrix(SA7.coords.rh), nrot = 10000) #rotate the SA7 parcellation 10,000 times on the freesurfer sphere to generate spatial nulls for spin-based permutation significance testing 
saveRDS(perm.id.full, "/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/SA7_sphericalrotations_N10000.rds")

