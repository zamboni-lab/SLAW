#THis part needs to be replaced by imports
#' @include distances.R
#' @include classes.R
#' @include allGenerics.R
#' @include utils.R

initializeParameters <- function(lar){
  lar@parameters <- list(rt=0.1,ppm=12,dmz=0.007,int=1,col_mz=1,rt_dens=0.04,
                         col_rt=2,col_int=3,sep_table=",",alignment_threshold=8, aligner_type="density")
  lar
}

#' Set the parameters on LCMSAlignReferences object
#'
#' @param lar An LCMSAlignerReferencesobject
#' @return A named list of the parameters of the dataset
#' @export
#'
#' @examples
#' print("Examples to be put here")
setMethod("params","LCMSAlignerReferences",function(lar){
  return(lar@parameters)
})


#' Set the parameters on LCMSAlignReferences object
#'
#' @param lar An LCMSAlignerReferencesobject
#' @param values A named list with names corresponding to the parameters.
#' @return A filled LCMSAlignerReference sobject
#' @export
#'
#' @examples
#' print("Examples to be put here")
setMethod("params<-","LCMSAlignerReferences",function(lar,values){
  pm <- match(names(values),names(lar@parameters))
  pna <- which(is.na(pm))
  if(length(pna)!=0){
    stop("Parameters ",paste(names(values)[pna],collapse=",")," do not exist.")
  }else{
    for(np in names(values)){
      lar@parameters[[np]] <- values[[np]]
    }
  }
  return(lar)
})





###eveentaully pick the same pick in every files.
LCMSAlignerReferences <- function(peaktables,num_peaks=150,col_int=character(0),
                                  col_mz=character(0),col_rt=character(0),full_interval = TRUE,
                                  bpp=NULL,seed=512,graphical=graphical,...){
  lar <- new("LCMSAlignerReferences")
  lar <- initializeParameters(lar)


  if(!is.null(bpparam)){
    lar@bpp <- bpp
  }else{
    lar@bpp <- bpparam()
  }

  ####We find the correct columns names
  phead <- getHeader(peaktables[[1]],intname = col_int,rtname = col_rt,mzname = col_mz)
  cpar <- as.list(phead)
  names(cpar) <-  c("col_mz","col_rt","col_int")
  params(lar) <- cpar

  ###Now we remove the rest
  ###We transform all the peaktables.
  peaktables <- sapply(peaktables,fortifyPeaktable,lar=lar,simplify=FALSE)


  eps <- extractWellBehavedPeaksManyFiles(lar,peaktables, sel_peaks_max=num_peaks,full_interval=full_interval,
                                          seed=seed,bpp=bpp,graphical=graphical,...)

  p_int <- match(col_int,colnames(eps[[1]]))

  vpeaks <- sapply(eps,initializeModelPeak,simplify=FALSE,rt_dev=lar@parameters$rt)
  temp <- do.call(rbind,vpeaks)

  idx_samples <- matrix(0,nrow=nrow(temp),ncol=length(peaktables))
  bins <- c(0,cumsum(sapply(peaktables,nrow)))
  bins[length(bins)] <- bins[length(bins)]+1

  vbins <- sapply(eps,function(x,bins){
    .bincode(as.numeric(row.names(x)),breaks = bins)
  },simplify=FALSE,bins=bins)

  for(i in seq_along(vbins)){
    idx_samples[i,vbins[[i]]] <- as.numeric(row.names(eps[[i]]))-bins[vbins[[i]]]
  }

  cfiles <- paste("idx_file_",1:length(peaktables),sep="")
  temp <- cbind(temp,idx_samples)
  ###We add the peak coordinates in each samples
  colnames(temp) <- c("mz","rt",col_int,"dmz","drt","dint","num",cfiles)
  lar@peaks <- data.frame(temp)
  message("Found ",nrow(lar@peaks)," references peaks.")
  return(lar)
}

initializeModelPeak <- function(peaks,ppm = 10, dmz= 0.003, rt_dev = 10,int_dev = 0.3,factor=3){
  vmean <- apply(peaks,2,mean)
  tolmz <- mztol(mz = vmean[1],ppm = ppm, dmz = dmz)
  tolrt <- rt_dev
  tolint <- vmean[3]*0.5
  return(c(vmean,tolmz,tolrt,tolint,nrow(peaks)))
}


createDistance <- function(lar){
  to_return <- function(p1,p2,ppm=lar@parameters$ppm,rt=lar@parameters$rt,int=lar@parameters$int){
    divergence(p1,p2,ppm=lar@parameters$ppm,rt=lar@parameters$rt,int=lar@parameters$int)
  }
  return(to_return)
}

##This return a list of peaks, each of them with a nieghbours and time complexity

###We defined a scoring systme fo revery metabolites.
###n N*(0.01-diff(range(mzs))+)
extractWellBehavedPeaksManyFiles <- function(lar,peaktables,
                                             int_quantile=0.6, sel_peaks_max=150,full_interval = TRUE,
                                             points_by_files = 1500, n_clusters = 10, points_by_cluster=7,seed=512,
                                             graphical=FALSE,bpp=NULL){
  fdist <- createDistance(lar)

  if(is.null(bpp)) bpp <- bpparam()
  n_clusters <- (n_clusters-2)

  ###We just peak the int_qunatiles most intense peaks in every files
  peaktables <- lapply(peaktables,function(x,prob){
    vquant <- quantile(x[,3],probs=int_quantile)
    x[x[,3]>vquant,,drop=FALSE]
  },prob=int_quantile)

  ###We put every peaks of every samples in an Rtree
  seqsize <- sapply(peaktables,nrow)
  btab <- do.call(rbind,peaktables)
  vintbtab <- btab[,3]


  tab_clustering <- sapply(peaktables,function(x,num,int_quantile){
    oint <- order(x[,3],decreasing=TRUE)
    x[oint[1:min(nrow(x),num)],1:2]
  },num=points_by_files,int_quantile=int_quantile,simplify = FALSE)

  tab_clustering <- do.call(rbind,tab_clustering)
  rtr <- range(tab_clustering[,2])

  ###First 5% is always put in a single cluster
  rt_spe <- c((rtr[1]+diff(rtr)*0.05),(rtr[2]-diff(rtr)*0.05))

  tab_clustering <- tab_clustering[(tab_clustering[,2]>rt_spe[1])&
                                     (tab_clustering[,2]<rt_spe[2]),,drop=FALSE]
  
  if(nrow(tab_clustering)>n_clusters){

  meanv <- apply(tab_clustering,2,mean)
  sdv <- apply(tab_clustering,2,sd)

  tab_clustering[,1] <- (tab_clustering[,1]-meanv[1])/sdv[1]
  tab_clustering[,2] <- (tab_clustering[,2]-meanv[2])/sdv[2]
  set.seed(seed)
  gmm <- GMM(tab_clustering,gaussian_comps = n_clusters,seed=seed)

  rt_range <- range(btab[,2])

  dtf <- data.frame(mz=tab_clustering[,1]*sdv[1]+meanv[1],rt=tab_clustering[,2]*sdv[2]+meanv[2],
                    cluster = apply(gmm$Log_likelihood,1,which.max))

    if(graphical)  plot(dtf$rt,dtf$mz,ylab="m/z",xlab="Retention time",col=as.factor(dtf$cluster),pch=16)
  
  }

  row.names(btab) <- NULL
  bfiles <-   c(0,cumsum(sapply(peaktables,nrow)))
  rtree <- RTree(as.matrix(btab[,c(1,2)]))
  k <- as.integer(length(peaktables))

  ####We first keep the most intense peaks only
  candidates <- which((vintbtab>quantile(vintbtab,1-int_quantile)))
  tcandidates <- as.matrix(btab[candidates,])

  ###We find the neughbors of the candidate points and keep only those present in each file
  neighbors <- knn(rtree,y = as.matrix(tcandidates[,1:2]),k=k)
  numfiles <- sapply(neighbors,.bincode,breaks=bfiles,include.lowest=TRUE,simplify=FALSE)
  num_present <- sapply(numfiles,function(x){length(unique(x))})
  opp <- order(num_present,decreasing=TRUE)

  # message("opp size: ",length(opp))

  # always_present <- which(sapply(numfiles,function(x,num,j){(num-length(unique(x)))<=j},num=length(peaktables),
  #                                j = allowed_jump))
  #
  # ###We keep each cluster a single time
  # neighbors_hash <- sapply(neighbors[always_present],function(x){
  #   paste(sort(x),collapse = "|")
  # })
  neighbors_hash <- sapply(neighbors,function(x){
    paste(sort(x),collapse = "|")
  })

  ###WE sort the cnadidates

  thash <- table(neighbors_hash)
  ohash <- order(thash,decreasing=TRUE)[1:min(length(thash),5000)]
  ###To avoid the computationall overlaod we only keep the 5000 first candidates


  #evident <- names(thash)[which(thash>=(length(peaktables)-allowed_jump))]
  ##Avoigin the raw filtration
  evident <- names(thash)[ohash]
  new_clusters <- match(evident,neighbors_hash)
  ####remaining candidates after this step.
  # sel_candidates <- always_present[new_clusters]
  candidates <- candidates[new_clusters]
  neighbors <- neighbors[new_clusters]

  # message("candidates size: ",length(neighbors))

  ####Now we refine the candidates by checking in each file
  rneighbors <- unlist(neighbors)
  class_check <- .bincode(rneighbors,bfiles)

  tocalculate <- split(rneighbors,drop=FALSE,f=class_check)

  ###we correct for the cumulative sum

  ####We calculate the distance to the most neighbouring point.
  inter_distances <- mapply(peaktables,tocalculate,bfiles[1:(length(bfiles)-1)],FUN=function(pt,ct,correction,k,
                                                                                             fdist){
    library(rtree)
    ct <- ct-correction
    rtr <- RTree(as.matrix(pt[,1:2]))
    tpt <- pt[ct,,drop=FALSE]
    kneigh <- knn(rtr,as.matrix(tpt[,1:2,drop=FALSE]),k=as.integer(length(peaktables)))

    ##We remove the peak himself
    kneigh <- mapply(kneigh,ct,FUN=function(x,y){
      setdiff(x,y)
    },SIMPLIFY = FALSE)

    vmap <- onlineLCMSaligner:::calcDistance(pt,tpt,kneigh,vfun=fdist)

    # vmap <- mapply(split(tpt,1:nrow(tpt)),kneigh,FUN=calcDistance,MoreArgs=list(peaktable=pt),SIMPLIFY = FALSE)
    return(vmap)
  },MoreArgs = list(k=k+1,fdist=createDistance(lar)),SIMPLIFY = FALSE)

  ####We  reformat all the data to allow a processing despit some missing samples
  idxv <- split(1:length(class_check),f = as.factor(class_check))
  reformatted_inter_distances <- rep(0,length(class_check))

  for(i in 1:length(idxv)){
    reformatted_inter_distances[idxv[[i]]] <- inter_distances[[i]]
  }

  reformatted_inter_distances <- split(as.data.frame(
    matrix(reformatted_inter_distances,ncol=k,byrow=TRUE)),f =
      1:(length(reformatted_inter_distances)/k))

  ###we aggregate the distance by min
  reformatted_inter_distances <- sapply(reformatted_inter_distances,min)


  ###We compute the minmmum distance of each cluster
  max_intra_dist <- lapply(neighbors,function(x,peaktable,fdist){
    # sel <- split(peaktable[x,,drop=FALSE],1:length(x))

    to_calc <- combn(x,2)

    vdist <- apply(to_calc,2,function(x,pt,fdist){
      fdist(pt[x[1],],pt[x[2],])
    },pt=peaktable,fdist=fdist)

    max(vdist)

  },peaktable=btab[,c(1,2,3)],fdist=fdist)
  max_intra_dist <- do.call(c,max_intra_dist)
  vratio <- as.numeric(reformatted_inter_distances)/as.numeric(max_intra_dist)
  o_ratio <- order(vratio,decreasing=TRUE)
  # message("oratio",length(o_ratio))
  sel_peaks <- sapply(neighbors,function(x,pt){pt[x,,drop=FALSE]},pt=btab,simplify = FALSE)

  ###We filter out all the peak which does not correspond to real group.
  mz_accurate <- sapply(sel_peaks,function(x,dmz,ppm){
    diff(range(x[,1]))/mean(x[,1])
  },dmz=lar@parameters$dmz,ppm=lar@parameters$ppm)
  # sel_peaks <- sel_peaks[mz_accurate]

  omz <- order(mz_accurate,decreasing=TRUE)

  ###We reorder the cluster based on the normalized ranking
  final_order <- order(scale(omz)[,1]+2*scale(ohash)[,1]+scale(o_ratio)[,1],decreasing=FALSE)
  sel_peaks <- sel_peaks[final_order]

  ###Likehood to a constant factor
  calc_likehood <- function(x,mn,sdn,cov_mat,centroids){
    x <- x[,c(1,2)]
    for(i in 1:ncol(x)){
      x[,i] <- (x[,i]-mn[i])/sdn[i]
    }

    like <- sapply(1:nrow(cov_mat),function(cc,peaks,cov_mat,centroids){
      covv <- cov_mat[cc,]
      centroids <- centroids[cc,]
      ###1st version the distance from the cluster centroid and that it
      vall <- apply(peaks,1,function(pp,covv,centroid){
        -0.5*(log(sum(covv))+matrix(pp-centroid,ncol=2)%*% diag(1/covv) %*% (pp-centroid))
      },covv=covv,centroid=centroids)
      sum(vall)
    },cov_mat=cov_mat,centroids=centroids,peaks=x)
    return(like)
  }
  
  mean_rt <-  sapply(sel_peaks,function(x){mean(x[,2])})
  clusters <- rep(0,length(sel_peaks))
  if(nrow(tab_clustering)>n_clusters){
  clusters <- sapply(sel_peaks,calc_likehood,mn=meanv,sdn=sdv,
         cov_mat=gmm$covariance_matrices,centroids=gmm$centroids,simplify=FALSE)
  clusters <- sapply(clusters,which.max)
  ##Low intensities points are put in theri own cluster
  }else{
    n_clusters <- 0
  }
  
  clusters[mean_rt<rt_spe[1]] <- n_clusters+1
  clusters[mean_rt>rt_spe[2]] <- n_clusters+2

  ###We keep at leat k points eventually
  points_by_cluster <- ceiling(sel_peaks_max/(n_clusters+2))

  count_bins <- rep(0,n_clusters+2)
  sel_idx <- rep(0,(n_clusters+2)*points_by_cluster)
  clustv <- rep(0,(n_clusters+2)*points_by_cluster)
  num_sel <- 0
  for(i in seq_along(clusters)){
    if(count_bins[clusters[i]]<=points_by_cluster){
      count_bins[clusters[i]] <- count_bins[clusters[i]]+1
      num_sel <-  num_sel + 1
      sel_idx[num_sel] <- i
      clustv[num_sel] <- clusters[i]
    }
  }

  ###We do a second pass to complete the data with the most probably ones.
  for(i in seq_along(clusters)){
    if(num_sel>=sel_peaks_max){
      break
    }
    if(!(i %in% sel_idx)){
      num_sel <-  num_sel + 1
      sel_idx[num_sel] <- i
    }
  }
  sel_idx <- sel_idx[1:num_sel]
  clustv <- clustv[1:num_sel]
  sel_peaks <- sel_peaks[sel_idx]

  return(sel_peaks)
}



# options(scipen = 0.0001)

setMethod("plot","LCMSAlignerReferences",function(x,y=NULL){
  df <- 2
})

getBoxLimit <- function(lar,model_peak){
  dx <- mztol(mz = model_peak[1],dmz = lar@parameters$dmz,ppm = lar@parameters$ppm)
  dy <- ifelse(length(model_peak)>=4,model_peak[4],1)*lar@parameters$rt
  return(c(dx,dy))
}
