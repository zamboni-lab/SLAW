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
      inboxes[[i]] <- withinBox(rtree_data,npoints[i,c(1,2),drop=FALSE],dx=dev_mz[i],dy=3*dev_rt)[[1]]+1
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


###the RT is estimated as a RANSAC LOESS AS We don t want a mistake
estimatesRTDeviationModel <- function(lar,peaktable,rt_scaling = c(0.9,0.95,1,1.05,1.1),
                                      ransac_niter=2000,ransac_dist_threshold=0.1,ransac_l1=0.05,
                                      span=span,frac_found=0.3,graphical=FALSE){
  if(!is.matrix(peaktable)) peaktable <- as.matrix(peaktable[,c(1,2,3)])

  rtree <- RTree(peaktable[,c(1,2)])
  rmatch <- findReferencesPoints(lar,rtree,peaktable,rt_scaling = rt_scaling)

  ###We remove the unmatched peaks
  rloess <- NULL
  non_matched <- is.na(rmatch)
  psel <- which(!non_matched)

  if((length(psel)/nrow(lar@peaks))<frac_found){
    warning("No enough references points found, no correction is done")
    ###We don t correct anything
    x <- 1:10
    y <- rep(0,10)
    rloess <- lm(y ~ x)
  }else{
  deviations <- (lar@peaks[psel,c(1,2)]-peaktable[rmatch[psel],c(1,2)])

  ####A first correction is always applied as the number of peak increase
  ###We now return the deviation.
  df_ransac <- data.frame(x=peaktable[rmatch[psel],2],y=deviations[,2])

  ###RObuest regression
  num_model <- nrow(df_ransac)

  # rloess <- tryCatch(ransacLoess(df_ransac,max_iter = ransac_niter,fitting_point = ceiling(num_model/2),
  #                       min_inliers = ceiling(num_model*1/4),dist_threshold = ransac_dist_threshold,L1=ransac_l1,span=span,plotL=TRUE),
  #                    error=function(e){return(NA)})
  rloess <- ransacLoess(df_ransac,max_iter = ransac_niter,fitting_point = ceiling(num_model/2),
                                 min_inliers = ceiling(num_model*1/4),dist_threshold = ransac_dist_threshold,L1=ransac_l1,span=span,
                                 graphical=graphical)

  }
  # cdist <- ransac_dist_threshold
  # it <- 1
  # while(it<2&(is.null(rloess)||is.na(rloess))){
  #   it <- it+1
  #   cdist <- cdist * 2
  #   rloess <- tryCatch(ransacLoess(df_ransac,max_iter = ransac_niter,fitting_point = ceiling(num_model/3),
  #                                  min_inliers = ceiling(num_model*1/2),span=span,dist_threshold = cdist),
  #                      error=function(e){return(NA)})
  # }

  ###We dont ocrrect anything os linear model.
  if(is.null(rloess)||is.na(rloess)){
    x <- 1:10
    y <- rep(0,10)
    rloess <- lm(y ~ x)
  }


  if((length(rloess)==1) && is.na(rloess)){
    warning("LOESS model foes not fit well enough try to tune the ransac_dist_threshold parameter")
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

  # too_low_rts <- which(peaktable[,2]<(pred_rt_range[1]-inter_to_predict))

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

  if(any(abs(cpeaktable-peaktable[,2]))>lim) peaktable[,2] <- peaktable[,2]

  return(peaktable)
}

correctPeaktable <- function(lar,peaktable,rt_scaling = c(0.9,0.95,1,1.05,1.1),
                             ransac_niter=4000,ransac_dist_threshold=0.1,ransac_l1=0.05,extensions=c(4,0.15),
                             span=0.6,lim=0.1,graphical=FALSE){

  rtmodel <- estimatesRTDeviationModel(lar,peaktable,rt_scaling = rt_scaling,ransac_niter=ransac_niter,
                            ransac_dist_threshold=ransac_dist_threshold,span=span,ransac_l1=ransac_l1,
                            graphical=graphical)

  corrected_peaktable <- correctRt(rtmodel,peaktable,extensions=extensions,lim=lim)


  return(list(peaktable=corrected_peaktable,model = rtmodel))
}
