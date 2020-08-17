getConsensusPoint <- function(peaktable,points,neighbours){
  vmap <- mapply(split(points,1:nrow(points)),neighbours,FUN=function(x,y){
    vmz <- peaktable[y[which.min(abs(peaktable[y,1]-x[1]))],1]
    vrt <- peaktable[y[which.min(abs(peaktable[y,2]-x[2]))],2]
    vint <- peaktable[y[which.min(abs(peaktable[y,3]-x[3]))],3]
    return(c(vmz,vrt,vint))
  },SIMPLIFY = FALSE)
  return(vmap)
}

###
divergence <- function(p1,p2,ppm=5,rt=5,int=1){
  ###foruml is deviation 
  p2 <- as.numeric(p2)
  p1 <- as.numeric(p1)
  vdiff <- p2-p1
  ppm_diff <- vdiff[1]*min(p2[1],p1[1])/1e6
  rt_diff <- abs(vdiff[2])/rt
  int_diff <- 0
  if((length(p1)>=3)&(length(p2)>3)){
  int_diff <- abs(vdiff[3])/max(p1[3],p2[3])*int
  }
  res <- ppm_diff+rt_diff+int_diff
  return(res)
}



calcDivergenceNN <- function(ref,points,neighbours,...){
  apply(X = points[neighbours,,drop=FALSE],MARGIN = 1,FUN = divergence,p1=ref,...)
}



calcDistance <- function(peaktable,candidates,neighbours,vfun=divergence){
  cpoints <- getConsensusPoint(peaktable,candidates,neighbours)
  mapply(split(as.data.frame(candidates),
               1:nrow(candidates)),cpoints,FUN=vfun,SIMPLIFY = TRUE)
}



