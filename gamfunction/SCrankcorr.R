library(R.matlab)
library(psych)
source("/Users/xuxiaoyu_work/Cuilab/GeneralRfunctions/rotate_parcellation/perm.sphere.p.R")
source("/Users/xuxiaoyu_work/Cuilab/DMRI_network_development/SC_development/Rcode_SCdevelopment/gamfunction/permspheregam.R")

SCrankcorr <- function(gamresult, computevar, ds.resolution, perm.id.full, dsdata=FALSE,
                       setgam=FALSE, set_fx=FALSE){

  #### connectional rank
  ########################
  #connectional rank
  Matrix.ds<-matrix(NA, nrow=ds.resolution, ncol=ds.resolution)
  indexup.ds <- upper.tri(Matrix.ds)
  indexsave.ds <- !indexup.ds
  Matrix.ds.index<-Matrix.ds
  SClength.ds=ds.resolution*(ds.resolution+1)/2
  Matrix.ds.index[indexsave.ds]<-c(1:SClength.ds)
  Matrix.ds.SCrank<-Matrix.ds
  for (x in 1:ds.resolution){
    for (y in 1:ds.resolution){
      Matrix.ds.SCrank[x,y]<-x^2+y^2
    }
  }
  Matrix.ds.SCrank[indexup.ds]<-NA
  Matrix.ds.SCrank[indexsave.ds]<-rank(Matrix.ds.SCrank[indexsave.ds], ties.method = "average")
  gamresult.ds<-data.frame(SCrank=Matrix.ds.SCrank[indexsave.ds], computevar=NA)
  gamresult.ds$computevar <-gamresult[, computevar]
  
  correstimate<-corr.test(gamresult.ds, method="spearman")$r[2,1]
  
  ## spin test
  
  perm.id.ds.L <- perm.id.full[1:12, ]
  perm.id.ds.R <- perm.id.full[13:24, ]-12
  perm.id.ds <- cbind(perm.id.ds.L, perm.id.ds.R)
  
  SC.perm.id.ds<-matrix(data=NA, nrow = SClength.ds, ncol=20000)
  SC.perm.id.sep<-mclapply(1:20000, function(i){
    ordertmp<-perm.id.ds[,i]
    tmp<-Matrix.ds.index[ordertmp, ordertmp]
    for (x in 1:ds.resolution){
      for (y in 1:ds.resolution){
        if (is.na(tmp[x,y])){
          tmp[x,y]<-tmp[y,x]}
      }}
    return(tmp[indexsave.ds])}, mc.cores = 4)
  for (j in 1:20000){
    SC.perm.id.ds[,j]<-SC.perm.id.sep[[j]]
  }
  ## if NA
  if (length(which(is.na(gamresult.ds$computevar)))>0){
    SC.perm.id.ds <- SC.perm.id.ds[-which(is.na(gamresult.ds$computevar)),]
    SC.perm.id.ds.test<-lapply(c(1:20000), function(x) rank(SC.perm.id.ds[,x], ties.method="first"))
    for (i in 1:20000){
      SC.perm.id.ds[,i]<-SC.perm.id.ds.test[[i]]
    }
    gamresult.ds.noNA <- gamresult.ds[-which(is.na(gamresult.ds$computevar)), ]
  }else{gamresult.ds.noNA <- gamresult.ds}
  
  if (setgam==FALSE){
    pspin <- perm.sphere.p(gamresult.ds.noNA$computevar, gamresult.ds.noNA$SCrank, SC.perm.id.ds,corr.type='spearman')
  }else{
    pspin <- perm.sphere.p.gam(gamresult.ds.noNA$SCrank, gamresult.ds.noNA$computevar, SC.perm.id.ds, set_fx)
  }
  
  SCrankdf <- data.frame(ds.resolution=ds.resolution, Interest.var=computevar,
                         r.spearman=correstimate, p.spin=pspin)
  if (dsdata == TRUE){
    names(gamresult.ds) <- c("SCrank", computevar)
    return(gamresult.ds)
  }else{
    return(SCrankdf)
  }
  
}



