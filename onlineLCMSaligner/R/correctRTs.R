#' @include classes.R
#' @include allGenerics.R
#' @include utils.R
#' @include ransac.R

###We remove all the value of this one
findReferencesPoints <- function(lar,rtree_data,peaktable_data,rt_scaling = 1){
  points <- as.matrix(lar@peaks[,c(1,2,3)])
  fdist <- createDistance(lar)

  ####We find all the neighbouring points using mz and rt
  dev_mz <- sapply(points[,1],mztol,ppm=lar@parameters$ppm,dmz=lar@parameters$dmz)
  dev_rt <- lar@parameters$rt
  dev_int <- lar@parameters$int

  ###WE test multiple deviation and score them if needed.
  best_matching <- NULL
  best_error <- 1e9

  for(ldev in rt_scaling){
    npoints <- points
    npoints[,2] <- npoints[,2]*ldev
    npoints <- as.matrix(npoints)
    inboxes <- vector(mode="list",length=nrow(npoints))
    for(i in seq_along(inboxes)){
      inboxes[[i]] <- withinBox.RTree(rtree_data,npoints[i,c(1,2),drop=FALSE],dx=dev_mz[i],dy=3*dev_rt)[[1]]+1
    }
    penalties <- mapply(split(npoints,f = 1:nrow(npoints)),dev_mz,FUN=function(x,dx,dy,fdist){
      sxp <- c(x[1]+dx,x[2]+5*dy,x[3]*3)
      fdist(x,sxp)
    },MoreArgs = list(fdist=fdist,dy=dev_rt))
    pmatch <- mapply(split(npoints,f = 1:nrow(npoints)),inboxes,penalties,FUN=function(x,y,pen,fdist,peaktable){
      if(length(y)==0){
        to_return <- c(NA,pen)
        return(to_return)
      }
      all_dist <- apply(peaktable[y,,drop=FALSE],1,fdist,p1=x)
      pclosest <- which.min(all_dist)
      to_return <- c(y[pclosest],all_dist[pclosest])
      return(to_return)
    },MoreArgs=list(fdist=fdist,peaktable=peaktable_data),SIMPLIFY=TRUE)


    ####We get the usm of the deviation
    current_error <- sum(pmatch[2,])
    if(best_error>current_error){
      best_matching <- pmatch[1,]
      best_error <- current_error
    }else{
      break
    }
  }
  return(best_matching)
}



identity_fun <- function(){
  x <- 1:10
  y <- rep(0,10)
  rloess <- lm(y ~ x)
  return(rloess)
}

###the RT is estimated as a RANSAC LOESS AS We don t want a mistake
estimatesRTDeviationModel <- function(lar,peaktable,peaks,rt_scaling = c(0.9,0.95,1,1.05,1.1),
                                      ransac_niter=10000,ransac_dist_threshold=0.1,ransac_l1=0.05,
                                      span=span,frac_found=0.3,ratio=1.0,graphical=FALSE){
  if(!is.matrix(peaktable)) peaktable <- as.matrix(peaktable[,c(1,2,3)])
  rtree <- RTree(peaktable[,c(1,2)])
  rmatch <- findReferencesPoints(lar,rtree,peaktable,rt_scaling = rt_scaling)
  
  ###We remove the unmatched peaks
  rloess <- NULL
  non_matched <- is.na(rmatch)
  psel <- which(!non_matched)
  
  maxDmz <- NULL
  maxDrt <- NULL
  matchedDrt <- NULL
  deviations <- NULL
  
  if((length(psel)/nrow(lar@peaks))<frac_found){
    warning("No enough references points found, default parameters used")
    ###We don t correct anything
    maxDmz <- lar@parameters$dmz
    maxDrt <- lar@parameters$rt*2
    matchedDrt <- max(lar@parameters$rt/2,diff(range(peaktable[,2]))/200)
  }else{
    deviations <- (lar@peaks[psel,c(1,2)]-peaktable[rmatch[psel],c(1,2)])
    
    ###We get the maximum deviation on references peaks 
    maxDmz <- max(median(abs(deviations[,1]))*3,0.002)
    maxDrt <- median(abs(deviations[,2]))*3
    matchedDrt <- max(median(abs(deviations[,2]))/2,diff(range(peaktable[,2]))/200)
  }
  
  # message("Init_par")
  ###We find the candidates
  df_ransac <- NULL
  num_points <- NULL
  # num_points <- ceiling(nrow(deviations)/2)
  min_inliers <- ceiling(num_points*1.3*ratio)
  unique_p <- 0
  if(nrow(peaks)>0){
    num_points <- ceiling(nrow(peaktable)/3)
    matched_mz <- xcms:::fastMatch(peaktable[,1],peaks[,1],tol=maxDmz)
    matched_rt <- xcms:::fastMatch(peaktable[,2],peaks[,2],tol=maxDrt)
    unique_p <- 0
    for(i in seq_along(matched_mz)){
      if(is.null(matched_mz[[i]])|is.null(matched_rt[[i]])) next
      tcand <- intersect(matched_mz[[i]],matched_rt[[i]])
      if(length(tcand)==0) next
      unique_p <- unique_p+1
      num_candidates <- length(tcand)
      matched_mz[[i]] <- cbind(rep(peaktable[i,2],num_candidates),peaks[tcand,2]-peaktable[i,2])
    }
    if(unique_p==nrow(peaktable)){
      return(identity_fun())
    }
    df_ransac <- do.call(rbind,matched_mz)
    colnames(df_ransac) <- c("x","y")
    df_ransac <- as.data.frame(df_ransac)
    num_points <- ceiling(nrow(peaktable)/3)
    min_inliers <- ceiling(min(floor(nrow(peaktable)/4),nrow(df_ransac)/2)*ratio)
  }else{
    ####A first correction is always applied as the number of peak increase
    ###We now return the deviation.
    if(is.null(deviations)) return(identity_fun())
    num_points <- length(unique(rmatch[psel]))
    df_ransac <- data.frame(x=peaktable[rmatch[psel],2],y=deviations[,2])
    unique_p <- nrow(df_ransac)
  }
  # message("Next")
  ###In every case we remove the datapoint moving by more than 20% by more than 20%
  threshold <- 0.2*diff(range(peaktable[,2]))
  df_ransac <- df_ransac[abs(df_ransac[,2])<threshold,]
  

  ###RObuest regression
  num_model <- length(psel)
  min_inliers <- ceiling((unique_p-length(psel))/2)
  
  if(num_points>nrow(df_ransac)){
    warning("Not enough references point found. No correction is applied.")
    rloess <- NA
  }else{
    rloess <- loRansacLoss(df_ransac,max_iter = ransac_niter,fitting_point = num_points,
                          min_inliers = min_inliers,dist_threshold = matchedDrt,L1=ransac_l1,span=span,
                          graphical=graphical)
  }
  # message("Loess done")
  ###We dont ocrrect anything os linear model.
  if(is.null(rloess)||is.na(rloess)){
    return(identity_fun())
  }


  if((length(rloess)==1) && is.na(rloess)){
    warning("LOESS modeldid  not fit well enough.")
  }

  ####We need to be sure that the loess alway include on correction
  return(rloess)
}

##Lim is the limit of the correction.
correctRt <- function(rtmodel,peaktable,extensions=c(4,0.05),maxCor = NULL,lim = 0.1){
  cpeaktable <- peaktable[,2]
  ort <- peaktable[,2]
  ###We predict the deviation for every RT flags
  pred_dev <- predict(rtmodel,newdata=data.frame(x=peaktable[,2]))

  pna <- is.na(pred_dev)
  pred_rt_range <- range(peaktable[!pna,2],na.rm = TRUE)
  total_range <- range(peaktable[,2])

  pnna <- which(!is.na(pred_dev))

  ###In case of liner model.
  if(length(pnna)==nrow(peaktable)){
    peaktable[pnna,2] <- peaktable[pnna,2]+pred_dev[pnna]
    return(peaktable)
  }

  ###For peak outside the rt range
  inter_to_predict <- diff(pred_rt_range)*extensions[2]
  x_low <- seq(pred_rt_range[1],pred_rt_range[1]+inter_to_predict,length=extensions[1])
  model_extensions_low <- predict(rtmodel,x_low)
  c_low <- coef(lm(model_extensions_low~x_low))
  low_rts <- which((peaktable[,2]<pred_rt_range[1])&
                     (peaktable[,2]>(pred_rt_range[1]-inter_to_predict)))

  if(length(low_rts)>0){
  cor_low_rts <- c_low[2]*peaktable[low_rts,2]+c_low[1]
  peaktable[low_rts,2] <- peaktable[low_rts,2]+cor_low_rts
  }



  x_high <- seq(pred_rt_range[2]-inter_to_predict,pred_rt_range[2],length=extensions[1])
  model_extensions_high <- predict(rtmodel,x_high)

  c_high <- coef(lm(model_extensions_high~x_high))

  high_rts <- which(peaktable[,2]>pred_rt_range[2]&
                      peaktable[,2]<pred_rt_range[2]+inter_to_predict)

  if(length(high_rts)>0){
    cor_high_rts <- c_high[2]*peaktable[high_rts,2]+c_high[1]
    peaktable[high_rts,2] <- peaktable[high_rts,2]+cor_high_rts
  }

  ###If the deviation is bigger than maxCor in variance we reduced it.
  if(!is.null(maxCor)){
    psel <- which(abs(pred_dev[pnna])>maxCor)
    if(length(psel)>0){
      pred_dev[pnna[psel]] <- maxCor*sign(pred_dev[pnna[psel]])
    }
  }

  ####We find the range end of the loess mode
  peaktable[pnna,2] <- peaktable[pnna,2]+pred_dev[pnna]

  if(any(abs(cpeaktable-peaktable[,2])>lim)) peaktable[,2] <- peaktable[,2]

  return(peaktable)
}

correctPeaktable <- function(lar,peaktable,peaks,rt_scaling = c(0.9,0.95,1,1.05,1.1),
                             ransac_niter=2000,ransac_dist_threshold=0.1,ransac_l1=0.05,extensions=c(4,0.15),
                             ransac_span=0.6,lim=0.1,ratio = 1.0,graphical=FALSE){

  rtmodel <- estimatesRTDeviationModel(lar,peaktable,peaks=peaks,rt_scaling = rt_scaling,ransac_niter=ransac_niter,
                                       span=ransac_span,ransac_dist_threshold=ransac_dist_threshold,
                                       ransac_l1=ransac_l1,ratio=ratio,graphical=graphical)

  corrected_peaktable <- correctRt(rtmodel,peaktable,extensions=extensions,lim=lim)


  return(list(peaktable=corrected_peaktable,model = rtmodel))
}
