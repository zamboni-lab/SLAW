#' @include utils.R
#' @include classes.R
#' @include allGenerics.R
#' @include LCMSAlignerFileHandler.R
#' @include LCMSAlignerDataStorage.R


###Insertion in the tree
getModelRtree <- function(lam) {
  rtree <- RTree(as.matrix(lam@peaks[, c(1, 2)]))
  return(rtree)
}

getPeaktableRtree <- function(peaktable, name_int = "int") {
  rtree <- RTree(as.matrix(peaktable[, c(1, 2)]))
  return(rtree)
}

#' Find cndidates to put into a peaktable.
#'
#' @param lar The reference object of an LCMSAlignerModel
#' @param rtree The Rtree of LC-MS  map data.
#' @param peak The peak to be inserted.
#'
#' @return The indices of the candidates
#'
#' @examples
#' print("Examples to be put here")
findCandidates <- function(lar, rtree, peak) {
  boxlims <- getBoxLimit(lar, peak)
  inbox <-
    withinBox.RTree(rtree, matrix(peak[c(1, 2)],ncol=2), boxlims[1], boxlims[2])
  return(inbox)
}

#' Compute all the divergence values of a set of graph
#'
#' @param lar The reference object of an LCMSAlignerModel
#' @param neighbours The index of the neighbouring points
#' @param points A three column with mz,rt,int
#' @param ref_point The reference point in whichi to comapre this shit
#'
#' @return
#' @export
#'
#' @examples
computeDivergence <- function(lar, neighbours, points, ref_point) {
  cdiv <-
    calcDivergenceNN(
      points = points,
      neighbours = neighbours,
      ref = ref_point,
      rt = lar@parameters$rt,
      ppm = lar@parameters$ppm,
      int = lar@parameters$int
    )
  return(cdiv)
}


#' @export
setMethod("show","LCMSAlignerModel",function(object){
  cat("An LCMSAlignerModel object including ",nrow(object@peaks),"unique peaks in ",numFiles(object),".")
})


numFiles <- function(lam){
  length(lam@files)
}


####Return the idx of the neighbors if the point is found without that return NA

#' Fidning the best candidates amongs a set of neighbours
#'
#' @param lar The reference object of an LCMSAlignerModel
#' @param neighbours The index of the neighbouring points
#' @param points The full set of points the data table
#' @param ref_point Te reference point to consider
#' @param threshold The threshold in divergence distance to csondier
#'
#' @return
#'
#' @examples
findBestCandidate <-
  function(lar,
           neighbours,
           points,
           ref_point,
           threshold = 5) {
    div <- computeDivergence(lar, neighbours, points, ref_point)
    pmin <- which.min(div)
    if (div[pmin] > threshold) {
      return(NA)
    }
    return(neighbours[pmin])
  }





updateFeatureInfos <- function(lam, peaks, pos) {
  vp <- by(cbind(peaks,pos),INDICES = pos,FUN = function(x){
    c(sum(x[,1]*x[,3])/sum(x[,3]),sum(x[,2]*x[,3])/sum(x[,3]),nrow(x),x[1,ncol(x)])
  })
  vp <- do.call(rbind,vp)

  lam@peaks[vp[,4], "num"] <- (lam@peaks[vp[,4], "num"] + vp[,3])
  lam@peaks[vp[,4], "current_num"] <- (lam@peaks[vp[,4], "current_num"] + vp[,3])

  lam@peaks[vp[,4], 1] <-
    (lam@peaks[vp[,4], 1] * (lam@peaks[vp[,4], "num"] - vp[,3]) + vp[,1]*vp[,3]) /
    (lam@peaks[vp[,4], "num"])

  lam@peaks[vp[,4], 2] <-
    (lam@peaks[vp[,4], 2] * (lam@peaks[vp[,4], "num"] - vp[,3]) + vp[,2]*vp[,3]) /
    (lam@peaks[vp[,4], "num"])
  return(lam)
}



#'
#' #' Insert a feature into the information table.
#' #'
#' #' @param lam An LCMS aligner model
#' #' @param peak A peak information
#' #' @param pos THe position in whichi to sinsert the data
#' #'
#' #' @return The LCMSAlignerModel
#' #' @export
#' #'
#' #' @examples
#' insertPointInFeature <- function(lam, peak, cor_peak, pos, id_sample) {
#'   lam@data[[pos]] <- extendDfRow(lam@data[[pos]], lam@peaks[pos, "current_num"])
#'   lam@data[[pos]][lam@peaks[pos, "current_num"], ] <- as.numeric(c(peak[[1]], peak[[2]],cor_peak[[2]],peak[[3]],id_sample))
#'   return(lam)
#' }

#' Align a peaktable to a model. The index are returned because the majroity of the data are not considered enough
#'
#' @param lam An LCMSAlignerModel
#' @param peaktable A peak table to be added to the model
#' @param id_sample The name of the ID to use eventually
#' @param threshold The threshold in diveergence (5 should be fine :) used ot threshold the model
#'
#' @return Am LCMSAlignerModel with the added peaktable.
#' @export
#'
#' @examples
alignToModel.nn <- function(lam,
                         peaktable,id_sample,
                         threshold = 5) {
  peaktable <- fortifyPeaktable(peaktable,lam@references)
  if (!is.matrix(peaktable)) {
    peaktable <- as.matrix(peaktable)
  }

  ####Peak table intensity ordering
  om <- order(peaktable[, 3], decreasing = TRUE)
  rtref <- getModelRtree(lam)
  new_group <- rep(TRUE, nrow(peaktable))
  new_index <- rep(NA_real_, nrow(peaktable))
  is_used_data <- rep(FALSE,nrow(lam@peaks))
  ###Rtree data structure
  # rtree <- RTree(as.matrix(peaktable[, c(1, 2)]))

  ###we align the points to model
  cpercent <- 0
  message("Beginning alignment. ", appendLF = FALSE)

  ###TODO parallelize
  all_cands <- apply(peaktable,1,function(x,rtree,lref,pt){
    sel_peak <- x
    ###We find the candidates
    fc <- findCandidates(lref, rtree = rtree, peak = x)[[1]]
        return(fc)
  },lref = lam@references,rtree = rtref,pt=lam@peaks)

  duplicates <- seq_along(all_cands)
  nits <- 0
  repeat{
    nits <- nits + 1

    ###1 Find all the best candidates for eveyr metbaolites
    new_index[duplicates] <- sapply(duplicates,function(idx,candidates,model_peaks,peaktable,lref,threshold=threshold){
      if(length(candidates[[idx]])==0) return(NA)
      new_index <- findBestCandidate(lref,
                        candidates[[idx]],
                        model_peaks,
                        ref_point = peaktable[idx,],
                        threshold = threshold)

      return(new_index)
    },candidates=all_cands,model_peaks = lam@peaks[,c(1,2)],
    peaktable = peaktable,lref=lam@references, threshold = threshold)

    ###We map findDuplicatePosition
    npfeat <- findDuplicatePosition(new_index[duplicates])

    ###We get the duplicated elements
    pdup <- which(npfeat!=-1)

    if(length(pdup)==0) break

    for(i in pdup){
      all_cands[[duplicates[i]]] <- setdiff(npfeat[i],all_cands[[duplicates[i]]])
    }

    duplicates <- duplicates[pdup]
  }

  message("Alignement finished after ",nits," iterations.")
  return(new_index)
}


###Align the data using an online estimatio of the density
#TODO change the binsize to a more rational qunaitity depending of the resolution
alignToModel.density <- function(lam, peaktable,bw = 10,binSize = 0.01, maxFeatures = 50,
                                 sleep = 0, mzCost = 0.005, rtCost = 0.05, bpp = NULL){

  ref_peaks <- lam@peaks[lam@order_peaks,c("mz","rt","num"),drop=FALSE]

  max_len <- nrow(peaktable)
  peaktable <- fortifyPeaktableIndex(peaktable,lam@references)
  vindex <- peaktable[[2]]
  peaktable <- peaktable[[1]]
  if (!is.matrix(peaktable)) {
    peaktable <- as.matrix(peaktable)
  }

  if (missing(peaktable))
    stop("Parameter 'peaktable' is missing!")
  if (!is.matrix(peaktable) & !is.data.frame(peaktable)){
    stop("'peaks' has to be a 'matrix' or a 'data.frame'!")
  }
  groupindex <- rep(NA_real_,max_len)

  cost_function <- function(x,y,mzCost,rtCost){
    sum(abs(x-y)/c(mzCost,rtCost))
  }
  cost_function <- Vectorize(cost_function)

  colt <- c("time","rt")
  vrt <- match(colt,colnames(peaktable))
  vrt <- colnames(peaktable)[vrt[!sapply(vrt,is.na)][1]]
  .reqCols <- c("mz", vrt)
  if (!all(.reqCols %in% colnames(peaktable)))
    stop("Required columns ", paste0("'", .reqCols[!.reqCols %in%
                                                     colnames(peaktable)], "'", collapse = ", "), " not found in 'peaktable' parameter")

  # peakOrder <- order(peaktable[, "mz"])
  rownames(peaktable) <- NULL
  rtRange <- range(peaktable[, vrt])
  mass <- seq(peaktable[1, "mz"], peaktable[nrow(peaktable), "mz"] + binSize,
              by = binSize/2)

  ####We find the greatest higher integer.
  masspos <- xcms:::findEqualGreaterM(peaktable[, "mz"], mass)
  masspos[masspos>nrow(peaktable)] <- 0

  massposmod <- xcms:::findEqualGreaterM(ref_peaks[, 1], mass)
  massposmod[massposmod>nrow(lam@peaks)] <- 0

  ##We select the interesting to facillitate parallel computations.
  sel_idx <- which((masspos[3:length(mass)]-masspos[1:(length(mass)-2)])>0)

  ###Initialisation of the group matrix.
  densFrom <- rtRange[1] - 3 * bw
  densTo <- rtRange[2] + 3 * bw
  densN <- max(512, 2 * 2^(ceiling(log2(diff(rtRange)/(bw/2)))))
  endIdx <- 0
  num <- 0
  #
    
    sliceDensity <- function(startIdx,endIdx,startIdxMod,endIdxMod,peaktable,ref_peaks,
                                 bw,maxFeatures,vrt,densFrom,densTo,densN){
      if ((endIdx - startIdx) < 0)
        return(matrix(0,nrow=0,ncol=2))
      if ((endIdxMod - startIdxMod) < 0)
        return(matrix(0,nrow=0,ncol=2))
      
      curMat <- peaktable[startIdx:endIdx, , drop = FALSE]
      subModel <- ref_peaks[startIdxMod:endIdxMod, , drop = FALSE]
      
      ###We generate the density for the model
      den <- suppressWarnings(density(curMat[, vrt], bw = bw, from = densFrom,
                                      to = densTo, n = densN))
      
      ###We generate the density of the model
      denModel <- suppressWarnings(density(subModel[,2], bw = bw, from = densFrom,
                                           to = densTo, n = densN,weights = subModel[,3]))
      deny <- den$y*nrow(curMat)+sum(subModel[,3])*denModel$y
      maxy <- NULL
      maxden <- max(deny)
      snum <- 0
      sub_res <- matrix(0,nrow=0,ncol=2)
      vres <- vector(mode="list",length=maxFeatures)
      for(i in seq_along(vres)){
        vres[[i]] <- matrix(0,nrow=0,ncol=2)
      }
      num_feat <- 1
      while (deny[maxy <- which.max(deny)] > maxden/500 && snum <
             maxFeatures) {
        if(num_feat > maxFeatures) break
        grange <- xcms:::descendMin(deny, maxy)
        deny[grange[1]:grange[2]] <- 0
        peaktable_idx <- which(curMat[, vrt] >= den$x[grange[1]] &
                                 curMat[, vrt] <= den$x[grange[2]])
        
        ###We get the number of group which are matched together eventually.
        ref_idx <- which(subModel[, 2] >= den$x[grange[1]] &
                           subModel[, 2] <= den$x[grange[2]])
        snum <- snum + 1
        num <- num + 1
        
        ###Cas where on of the catheogry is missing under the density.
        if((length(peaktable_idx)==0)|(length(ref_idx)==0)) next
        #We update the value.
        selm <- sapply(curMat[peaktable_idx,2],function(x,ref){which.min(abs(ref-x))},
                       ref=subModel[ref_idx,2])
        
        
        vres[[num_feat]] <- matrix(c(startIdx+peaktable_idx-1,startIdxMod+ref_idx[selm]-1),
                                   ncol=2,byrow=FALSE)
        num_feat <- num_feat+1
      }
      
      if(num_feat>1){
        return(do.call(rbind,vres[1:num_feat]))
      }else{
        return(matrix(0,nrow=0,ncol=2))
      }
      
    }
    
    # COmputing the batch density by batch in differnet folder
    n <- bpworkers(bpp)
    if(is.list(n)){
      vapp <- mapply((masspos[1:(length(mass)-2)])[sel_idx],
                       (masspos[3:length(mass)])[sel_idx]-1,
                       (massposmod[1:(length(massposmod)-2)])[sel_idx],
                       (massposmod[3:length(massposmod)])[sel_idx]-1,
                       FUN=sliceDensity,MoreArgs=list(ref_peaks=ref_peaks,peaktable=peaktable,
                                                      bw=bw,maxFeatures=maxFeatures,vrt=vrt,
                                                      densFrom=densFrom,densTo=densTo,densN),SIMPLIFY = FALSE)
      
    }else{
      #We split in as many workers as possible
      startIdxs <- (masspos[1:(length(mass)-2)])[sel_idx]
      startIdxs <- split(startIdxs, rep_len(1:n, length(startIdxs)))
      endIdxs <- (masspos[3:length(mass)])[sel_idx]-1
      endIdxs <- split(endIdxs, rep_len(1:n, length(endIdxs)))
      startIdxMods <- (massposmod[1:(length(massposmod)-2)])[sel_idx]
      startIdxMods <- split(startIdxMods, rep_len(1:n, length(startIdxMods)))
      endIdxMods <- (massposmod[3:length(massposmod)])[sel_idx]-1
      endIdxMods <- split(endIdxMods, rep_len(1:n, length(endIdxMods)))
      sliceDensityGrouped <- function(startIdxs,endIdxs,startIdxMods,endIdxMods,supp_args){
        res <- mapply(startIdxs,endIdxs,startIdxMods,endIdxMods,
               FUN=sliceDensity,MoreArgs=supp_args,SIMPLIFY = FALSE)
        return(do.call(rbind,res))
      }
      vapp <- bpmapply(startIdxs,
                       endIdxs,
                       startIdxMods,
                       endIdxMods,
                       FUN=sliceDensityGrouped,MoreArgs=list(supp_args=list(ref_peaks=ref_peaks,peaktable=peaktable,
                                                                            bw=bw,maxFeatures=maxFeatures,vrt=vrt,
                                                                            densFrom=densFrom,densTo=densTo,densN)),
                       BPPARAM = bpp,SIMPLIFY = FALSE)
      vapp <- do.call(rbind,vapp)
      
    }
    

    # vapp <- bpmapply((masspos[1:(length(mass)-2)])[sel_idx],
    #                  (masspos[3:length(mass)])[sel_idx]-1,
    #                  (massposmod[1:(length(massposmod)-2)])[sel_idx],
    #                  (massposmod[3:length(massposmod)])[sel_idx]-1,
    #                  FUN=sliceDensity,MoreArgs=list(ref_peaks=ref_peaks,peaktable=peaktable,
    #                                  bw=bw,maxFeatures=maxFeatures,vrt=vrt,
    #                                  densFrom=densFrom,densTo=densTo,densN),
    #                  BPPARAM = bpp,SIMPLIFY = FALSE)
    
  
  ###parallel processing of the density estination to speed up the process.
  # vapp <- bpmapply((masspos[1:(length(mass)-2)])[sel_idx],
  #                  (masspos[3:length(mass)])[sel_idx]-1,
  #                  (massposmod[1:(length(massposmod)-2)])[sel_idx],
  #                  (massposmod[3:length(massposmod)])[sel_idx]-1,
  #          FUN=function(startIdx,endIdx,startIdxMod,endIdxMod,peaktable,ref_peaks,
  #                       bw,maxFeatures,vrt,densFrom,densTo,densN){
  #            if ((endIdx - startIdx) < 0)
  #              return(matrix(0,nrow=0,ncol=2))
  #            if ((endIdxMod - startIdxMod) < 0)
  #              return(matrix(0,nrow=0,ncol=2))
  # 
  #            curMat <- peaktable[startIdx:endIdx, , drop = FALSE]
  #            subModel <- ref_peaks[startIdxMod:endIdxMod, , drop = FALSE]
  # 
  #            ###We generate the density for the model
  #            den <- suppressWarnings(density(curMat[, vrt], bw = bw, from = densFrom,
  #                           to = densTo, n = densN))
  # 
  #            ###We generate the density of the model
  #            denModel <- suppressWarnings(density(subModel[,2], bw = bw, from = densFrom,
  #                                to = densTo, n = densN,weights = subModel[,3]))
  #            deny <- den$y*nrow(curMat)+sum(subModel[,3])*denModel$y
  #            maxy <- NULL
  #            maxden <- max(deny)
  #            snum <- 0
  #            sub_res <- matrix(0,nrow=0,ncol=2)
  #            vres <- vector(mode="list",length=maxFeatures)
  #            for(i in seq_along(vres)){
  #              vres[[i]] <- matrix(0,nrow=0,ncol=2)
  #            }
  #            num_feat <- 1
  #            while (deny[maxy <- which.max(deny)] > maxden/500 && snum <
  #                   maxFeatures) {
  #              if(num_feat > maxFeatures) break
  #              grange <- xcms:::descendMin(deny, maxy)
  #              deny[grange[1]:grange[2]] <- 0
  #              peaktable_idx <- which(curMat[, vrt] >= den$x[grange[1]] &
  #                                       curMat[, vrt] <= den$x[grange[2]])
  # 
  #              ###We get the number of group which are matched together eventually.
  #              ref_idx <- which(subModel[, 2] >= den$x[grange[1]] &
  #                                 subModel[, 2] <= den$x[grange[2]])
  #              snum <- snum + 1
  #              num <- num + 1
  # 
  #              ###Cas where on of the catheogry is missing under the density.
  #              if((length(peaktable_idx)==0)|(length(ref_idx)==0)) next
  #              #We update the value.
  #              selm <- sapply(curMat[peaktable_idx,2],function(x,ref){which.min(abs(ref-x))},
  #                             ref=subModel[ref_idx,2])
  # 
  # 
  #                vres[[num_feat]] <- matrix(c(startIdx+peaktable_idx-1,startIdxMod+ref_idx[selm]-1),
  #                                           ncol=2,byrow=FALSE)
  #                num_feat <- num_feat+1
  #            }
  # 
  #            if(num_feat>1){
  #              return(do.call(rbind,vres[1:num_feat]))
  #            }else{
  #              return(matrix(0,nrow=0,ncol=2))
  #            }
  # 
  # },MoreArgs=list(ref_peaks=ref_peaks,peaktable=peaktable,
  #                 bw=bw,maxFeatures=maxFeatures,vrt=vrt,
  #                 densFrom=densFrom,densTo=densTo,densN),
  # BPPARAM = bpp,SIMPLIFY = FALSE)


  # ####We only keep the best group
  # mapply(vapp[1:(length(vapp)-1)],vapp[1:(length(vapp)-1)],FUN = function(x,y){
  # tx <- table(x[,2])
  # ty <- table(y[,2])
  # ntx <- as.numeric(names(tx))
  # nty <- as.numeric(names(ty))
  # ###We check the intersection
  # inter <- intersect(x[,1],x[,2])
  # if(length(inter)!=0){
  #   vpos <- tx[as.character(inter)]>ty[as.character(inter)]
  #   ###We remove all the values eventually
  #   ux <- c(setdiff(ntx,inter),inter[!vpos])
  #   uy <- c(setdiff(ntx,inter),inter[vpos])
  # }
  #
  # },SIMPLIFY = FALSE,USE.NAMES = FALSE)
  #
  #
  # if (nrow(res)) {
  #   ## Remove groups that overlap with more "well-behaved" groups
  #   numsamp <- rowSums(
  #     as.matrix(res[, (match("npeaks", colnames(res)) +1):(ncol(res) -1),
  #                   drop = FALSE]))
  #   uorder <- order(-numsamp, res[, "npeaks"])
  #
  #   uindex <- rectUnique(
  #     as.matrix(res[, c("mzmin", "mzmax", "rtmin", "rtmax"),
  #                   drop = FALSE]), uorder)
  #   res <- res[uindex, , drop = FALSE]
  #   rownames(res) <- NULL
  # }
  #
  vsize <- lam@peaks[,"num"]
  vapp <- do.call(rbind,vapp)
  vapp <- bplapply(split(as.data.frame(vapp),f = vapp[,1],
                 drop = FALSE),FUN=function(x,vsize){x[which.max(vsize[x[,2]]),]},
           vsize=vsize,BPPARAM = bpp)
  vapp <- do.call(rbind,vapp)

  vtest <- tryCatch(nrow(vapp),error=function(e){
    return(NA)
  })
  if(nrow(vapp)==0) return(groupindex)

  ##If there is any NA we remove it (TODEBUG LATER)
  # vapp <- vapp[!is.na(vapp[,1]),,drop=FALSE]
  # vapp <- vapp[!duplicated(vapp[,1]),]
  groupindex[vindex[vapp[,1]]] <- vapp[,2]
  ###We just return the index
  return(groupindex)
}

alignToModelByBatch.density <- function(lam, peaktable,bw = 10,binSize = 0.01, maxFeatures = 50,
                                        sleep = 0, mzCost = 0.005, rtCost = 0.05, bpp = NULL){
  
  ref_peaks <- lam@peaks[lam@order_peaks,c("mz","rt","num"),drop=FALSE]
  
  max_index <- nrow(ref_peaks)+1
  
  opeaktable <- order(peaktable[,"mz"],decreasing=FALSE)
  
  peaktable <- peaktable[opeaktable,,drop=FALSE]
  
  groupindex <- rep(NA_real_,nrow(peaktable))
  cost_function <- function(x,y,mzCost,rtCost){
    sum(abs(x-y)/c(mzCost,rtCost))
  }
  cost_function <- Vectorize(cost_function)
  
  colt <- c("time","rt")
  vrt <- match(colt,colnames(peaktable))
  vrt <- colnames(peaktable)[vrt[!sapply(vrt,is.na)][1]]
  .reqCols <- c("mz", vrt)
  if (!all(.reqCols %in% colnames(peaktable)))
    stop("Required columns ", paste0("'", .reqCols[!.reqCols %in%
                                                     colnames(peaktable)], "'", collapse = ", "), " not found in 'peaktable' parameter")
  
  # peakOrder <- order(peaktable[, "mz"])
  rownames(peaktable) <- NULL
  rtRange <- range(peaktable[, vrt])
  mass <- seq(peaktable[1, "mz"], peaktable[nrow(peaktable), "mz"] + binSize,
              by = binSize/2)
  
  ####We find the greatest higher integer.
  masspos <- xcms:::findEqualGreaterM(peaktable[, "mz"], mass)
  masspos[masspos>nrow(peaktable)] <- 0
  
  massposmod <- xcms:::findEqualGreaterM(ref_peaks[, 1], mass)
  massposmod[massposmod>nrow(lam@peaks)] <- 0
  
  ##We select the interesting to facillitate parallel computations.
  sel_idx <- which((masspos[3:length(mass)]-masspos[1:(length(mass)-2)])>0)
  
  ##We increase the last bin to include the last point
  if(masspos[2+sel_idx[length(sel_idx)]]==nrow(peaktable)){
    masspos[2+sel_idx[length(sel_idx)]] <- nrow(peaktable)+1
  }
  
  ##We increase the last bin to include the last point
  if((massposmod[2+sel_idx[length(sel_idx)]]!=0) &
    (massposmod[2+sel_idx[length(sel_idx)]]==nrow(ref_peaks))){
    massposmod[2+sel_idx[length(sel_idx)]] <- nrow(ref_peaks)+1
  }
  
  ###Initialisation of the group matrix.
  densFrom <- rtRange[1] - 3 * bw
  densTo <- rtRange[2] + 3 * bw
  densN <- max(512, 2 * 2^(ceiling(log2(diff(rtRange)/(bw/2)))))
  endIdx <- 0
  num <- 0
  
  ###Slice processing function
  sliceDensity <- function(startIdx,endIdx,startIdxMod,endIdxMod,peaktable,ref_peaks,
           bw,maxFeatures,vrt,densFrom,densTo,densN,max_idx){
    if ((endIdx - startIdx) < 0)
      return(matrix(0,nrow=0,ncol=2))
    
    ##In this case ther can still be peak which are put in ther own groups
    ref <- TRUE
    if ((endIdxMod - startIdxMod) < 0)
      ref <- FALSE
    
    curMat <- peaktable[startIdx:endIdx, , drop = FALSE]
    
    if(ref) subModel <- ref_peaks[startIdxMod:endIdxMod, , drop = FALSE]
    
    ###We generate the density for the model
    den <- suppressWarnings(density(curMat[, vrt], bw = bw, from = densFrom,
                                    to = densTo, n = densN))
    
    ###We generate the density of the model
    if(ref){
      denModel <- suppressWarnings(density(subModel[,2], bw = bw, from = densFrom,
                                           to = densTo, n = densN,weights = subModel[,3]))
      deny <- den$y*nrow(curMat)+sum(subModel[,3])*denModel$y
    }else{
      deny <- den$y
    }
    maxy <- NULL
    maxden <- max(deny)
    snum <- 0
    sub_res <- matrix(0,nrow=0,ncol=2)
    vres <- vector(mode="list",length=maxFeatures)
    for(i in seq_along(vres)){
      vres[[i]] <- matrix(0,nrow=0,ncol=2)
    }
    num_feat <- 1
    while (deny[maxy <- which.max(deny)] > maxden/500 && snum <
           maxFeatures) {
      if(num_feat > maxFeatures) break
      grange <- xcms:::descendMin(deny, maxy)
      deny[grange[1]:grange[2]] <- 0
      peaktable_idx <- which(curMat[, vrt] >= den$x[grange[1]] &
                               curMat[, vrt] <= den$x[grange[2]])
      
      
      ###We get the number of group which are matched together eventually.
      if(ref){
        ref_idx <- which(subModel[, 2] >= den$x[grange[1]] &
                           subModel[, 2] <= den$x[grange[2]])
      }else{
        if(length(peaktable_idx)==1) next
      }
      snum <- snum + 1
      num <- num + 1
      
      ###Case where on of the catheogry is missing under the density.
      if(length(peaktable_idx)==0) next
      #We update the value.
      
      ###for each sample we select the closest retention time.
      if(ref && (length(ref_idx)>0)){
        selm <- sapply(curMat[peaktable_idx,2],function(x,ref){which.min(abs(ref-x))},
                       ref=subModel[ref_idx,2])
        
        
        vres[[num_feat]] <- matrix(c(startIdx+peaktable_idx-1,startIdxMod+ref_idx[selm]-1),
                                   ncol=2,byrow=FALSE)
      }else{
        vres[[num_feat]] <- matrix(c(startIdx+peaktable_idx-1,rep(max_idx+snum,length(peaktable_idx))),
                                   ncol=2,byrow=FALSE)
      }
      num_feat <- num_feat+1
    }
    
    if(num_feat>1){
      return(do.call(rbind,vres[1:num_feat]))
    }else{
      return(matrix(0,nrow=0,ncol=2))
    }
    
  }
  
  ###PArallel processing by split data.frame
  
  # COmputing the batch density by batch in differnet folder
  n <- bpworkers(bpp)
  message("Starting parallel processing ",Sys.time())
  if(is.list(n)){
    vapp <- mapply((masspos[1:(length(mass)-2)])[sel_idx],
                   (masspos[3:length(mass)])[sel_idx]-1,
                   (massposmod[1:(length(massposmod)-2)])[sel_idx],
                   (massposmod[3:length(massposmod)])[sel_idx]-1,
                   FUN=sliceDensity,MoreArgs=list(ref_peaks=ref_peaks,peaktable=peaktable,
                                                            bw=bw,maxFeatures=maxFeatures,vrt=vrt,
                                                            densFrom=densFrom,densTo=densTo,densN=densN,max_idx=max_index),SIMPLIFY = FALSE)
    
  }else{
    message("Grouped processing")
    #We split in as many workers as possible
    startIdxs <- (masspos[1:(length(mass)-2)])[sel_idx]
    startIdxs <- split(startIdxs, rep_len(1:n, length(startIdxs)))
    endIdxs <- (masspos[3:length(mass)])[sel_idx]-1
    endIdxs <- split(endIdxs, rep_len(1:n, length(endIdxs)))
    startIdxMods <- (massposmod[1:(length(massposmod)-2)])[sel_idx]
    startIdxMods <- split(startIdxMods, rep_len(1:n, length(startIdxMods)))
    endIdxMods <- (massposmod[3:length(massposmod)])[sel_idx]-1
    endIdxMods <- split(endIdxMods, rep_len(1:n, length(endIdxMods)))
    sliceDensityGrouped <- function(startIdxs,endIdxs,startIdxMods,endIdxMods,supp_args){
      res <- mapply(startIdxs,endIdxs,startIdxMods,endIdxMods,
                    FUN=sliceDensity,MoreArgs=supp_args,SIMPLIFY = FALSE)
      return(res)
    }
    vapp <- bpmapply(startIdxs,
                     endIdxs,
                     startIdxMods,
                     endIdxMods,
                     FUN=sliceDensityGrouped,MoreArgs=list(supp_args=list(ref_peaks=ref_peaks,peaktable=peaktable,
                                                                                    bw=bw,maxFeatures=maxFeatures,vrt=vrt,
                                                                                    densFrom=densFrom,densTo=densTo,densN=densN,max_idx=max_index)),
                     BPPARAM = bpp,SIMPLIFY = FALSE)
    vapp <- do.call(c,vapp)

    
  }
  message("End parallel processing ",Sys.time())
  
  
  # 
  # 
  # ###parallel processing of the density estination to speed up the process.
  # message("Starting parallel processing ",Sys.time())
  # vapp <- bpmapply((masspos[1:(length(mass)-2)])[sel_idx],
  #                  (masspos[3:length(mass)])[sel_idx]-1,
  #                  (massposmod[1:(length(massposmod)-2)])[sel_idx],
  #                  (massposmod[3:length(massposmod)])[sel_idx]-1,
  #                  FUN=function(startIdx,endIdx,startIdxMod,endIdxMod,peaktable,ref_peaks,
  #                               bw,maxFeatures,vrt,densFrom,densTo,densN,max_idx){
  #                    if ((endIdx - startIdx) < 0)
  #                      return(matrix(0,nrow=0,ncol=2))
  #                    
  #                    ##In this case ther can still be peak which are put in ther own groups
  #                    ref <- TRUE
  #                    if ((endIdxMod - startIdxMod) < 0)
  #                      ref <- FALSE
  #                    
  #                    curMat <- peaktable[startIdx:endIdx, , drop = FALSE]
  #                    
  #                    if(ref) subModel <- ref_peaks[startIdxMod:endIdxMod, , drop = FALSE]
  #                    
  #                    ###We generate the density for the model
  #                    den <- suppressWarnings(density(curMat[, vrt], bw = bw, from = densFrom,
  #                                   to = densTo, n = densN))
  #                    
  #                    ###We generate the density of the model
  #                    if(ref){
  #                      denModel <- suppressWarnings(density(subModel[,2], bw = bw, from = densFrom,
  #                                          to = densTo, n = densN,weights = subModel[,3]))
  #                      deny <- den$y*nrow(curMat)+sum(subModel[,3])*denModel$y
  #                    }else{
  #                      deny <- den$y
  #                    }
  #                    maxy <- NULL
  #                    maxden <- max(deny)
  #                    snum <- 0
  #                    sub_res <- matrix(0,nrow=0,ncol=2)
  #                    vres <- vector(mode="list",length=maxFeatures)
  #                    for(i in seq_along(vres)){
  #                      vres[[i]] <- matrix(0,nrow=0,ncol=2)
  #                    }
  #                    num_feat <- 1
  #                    while (deny[maxy <- which.max(deny)] > maxden/500 && snum <
  #                           maxFeatures) {
  #                      if(num_feat > maxFeatures) break
  #                      grange <- xcms:::descendMin(deny, maxy)
  #                      deny[grange[1]:grange[2]] <- 0
  #                      peaktable_idx <- which(curMat[, vrt] >= den$x[grange[1]] &
  #                                               curMat[, vrt] <= den$x[grange[2]])
  #                      
  #                      
  #                      ###We get the number of group which are matched together eventually.
  #                      if(ref){
  #                        ref_idx <- which(subModel[, 2] >= den$x[grange[1]] &
  #                                           subModel[, 2] <= den$x[grange[2]])
  #                      }else{
  #                        if(length(peaktable_idx)==1) next
  #                      }
  #                      snum <- snum + 1
  #                      num <- num + 1
  #                      
  #                      ###Case where on of the catheogry is missing under the density.
  #                      if(length(peaktable_idx)==0) next
  #                      #We update the value.
  #                      
  #                      ###for each sample we select the closest retention time.
  #                      if(ref && (length(ref_idx)>0)){
  #                        selm <- sapply(curMat[peaktable_idx,2],function(x,ref){which.min(abs(ref-x))},
  #                                       ref=subModel[ref_idx,2])
  #                        
  #                        
  #                        vres[[num_feat]] <- matrix(c(startIdx+peaktable_idx-1,startIdxMod+ref_idx[selm]-1),
  #                                                   ncol=2,byrow=FALSE)
  #                      }else{
  #                        vres[[num_feat]] <- matrix(c(startIdx+peaktable_idx-1,rep(max_idx+snum,length(peaktable_idx))),
  #                                                   ncol=2,byrow=FALSE)
  #                      }
  #                      num_feat <- num_feat+1
  #                    }
  #                    
  #                    if(num_feat>1){
  #                      return(do.call(rbind,vres[1:num_feat]))
  #                    }else{
  #                      return(matrix(0,nrow=0,ncol=2))
  #                    }
  #                    
  #                  },MoreArgs=list(ref_peaks=ref_peaks,peaktable=peaktable,
  #                                  bw=bw,maxFeatures=maxFeatures,vrt=vrt,
  #                                  densFrom=densFrom,densTo=densTo,densN=densN,max_idx=max_index),BPPARAM = bpp)
  # message("End parallel processing ",Sys.time())
  vsize <- lam@peaks[,"num"]
  ###if nothing has been aligned correctly.
  if(length(vapp)==0){
    return(groupindex)
  }
  
  ###We shift the index
  shift <- 0
  new_ids <- max_index+1
  temp_size <- rep(0,min(5000,2*nrow(ref_peaks)))
  if(nrow(ref_peaks)>0){
    temp_size[seq_along(vsize)] <- vsize
  }
  message("Starting refinement ",Sys.time())
  ### Vectorize version
  # if(FALSE){
    #tentative of vectorized processing of vapp
    
    # vapp <- list(
    #   matrix(c(1,2,4,2,2,5),ncol=2,byrow = FALSE),
    #   matrix(nrow=0,ncol=2),
    #   matrix(c(3,4,5,2,3,5),ncol=2,byrow=FALSE)
    # )
    # vsize <- c(2,1,4,3)
    # max_index <- 5
    ##We first have to count the number of shift
    limsize <- c(0,cumsum(sapply(vapp,nrow)))
    vapp <- do.call(rbind,vapp)
    vsupp <- which(vapp[,2]>=max_index)
    #We create a specific index
    #We then remove the rest
    modified_vsupp_bins <- .bincode(vsupp,limsize)
    vsupp_idx <- paste(modified_vsupp_bins,vapp[vsupp,2],sep="_")
    supp_ids <- unique(vsupp_idx)
    new_ids <- max_index+seq_along(supp_ids)
    temp_size <- numeric(new_ids[length(new_ids)])
    temp_size[seq_along(vsize)] <- vsize
    shifts <- rep(0,nrow(vapp))
    shifts[vsupp] <- match(vsupp_idx,supp_ids)
    vapp[,2] <- vapp[,2]+shifts
    vapp <- do.call(rbind,by(vapp,INDICES=vapp[,1],
                             FUN = function(x,vsize){x[which.max(vsize[x[,2]]),]},vsize=temp_size))
    psupp <- which(vapp[,2]>max_index)
    if(length(psupp)>=1){
      unique_values <- unique(vapp[psupp,2])
      vapp[psupp,2] <- match(vapp[psupp,2],unique_values)+max_index-1
      
    }
    
  # }
  
  ### LOOP VERSION
  # for(i in seq_along(vapp)){
  #   if(nrow(vapp[[i]])==0) next
  #   temp <- vapp[[i]]
  #   psupp <- which(temp[,2]>=max_index)
  #   ###Debuggued.
  #   pmodel <- which(temp[,2]<max_index)
  #   ##We add the size of the cluster to the new size.
  #   for(idx in temp[pmodel,2]){
  #     temp_size[idx] <- temp_size[idx]+1
  #   }
  #   if(length(psupp)==0){
  #     next
  #   }
  #   temp[psupp,2] <- temp[psupp,2]+shift
  #   
  #   ###We check the size of the new index
  #   if(max(temp[psupp,2])>length(temp_size)){
  #     temp_size <- c(temp_size,rep(0,length(temp_size)))
  #   }
  #   
  #   t2 <- table(temp[psupp,2])
  #   temp_size[temp[psupp,2]] <- t2[match(temp[psupp,2],names(t2))]
  #   vapp[[i]] <- temp
  #   shift <- shift+length(unique(temp[psupp,2]))
  # }
  # message("End refinement ",Sys.time())
  # vapp <- do.call(rbind,vapp)
  # vapp <- bplapply(split(as.data.frame(vapp),f = vapp[,1],
  #                        drop = FALSE),FUN=function(x,vsize){x[which.max(vsize[x[,2]]),]},
  #                  vsize=temp_size,BPPARAM = bpp)
  # 
  # vapp <- do.call(rbind,vapp)
  # ###We relabel with consecutive integer
  # psupp <- which(vapp[,2]>max_index)
  # if(length(psupp)>=1){
  #   unique_values <- unique(vapp[psupp,2])
  #   vapp[psupp,2] <- match(vapp[psupp,2],unique_values)+max_index-1
  #   
  # }
  message("End refinement ",Sys.time())
  if(nrow(vapp)==0) return(groupindex)
  ##If there is any NA we remove it (TODEBUG LATER)
  groupindex[opeaktable[vapp[,1]]] <- vapp[,2]
  ###We just return the index
  return(groupindex)
}



###A final clustering is odne to correct mistakes of the previous step eventually.
finalClustering <- function(lam,bw,binSize=0.01,bpp=NULL){

  if(is.null(bpp)) bpp <- bpparam()
  
  ###if a bug is constated we rresort the integer
  
  

  peaktable <- lam@peaks[lam@order_peaks,c("mz","rt","num"),drop=FALSE]


  # peakOrder <- order(peaktable[, "mz"])
  rownames(peaktable) <- NULL
  rtRange <- range(peaktable[, 2])
  mass <- seq(peaktable[1, "mz"], peaktable[nrow(peaktable), "mz"] + binSize,
              by = binSize/2)

  ####We find the greatest higher integer.
  masspos <- xcms:::findEqualGreaterM(peaktable[, "mz"], mass)
  masspos[masspos>nrow(peaktable)] <- 0

  ##We select the interesting to facillitate parallel computations.
  sel_idx <- which((masspos[3:length(masspos)]-masspos[1:(length(masspos)-2)])>0)

  ###Initialisation of the group matrix.
  densFrom <- rtRange[1] - 3 * bw
  densTo <- rtRange[2] + 3 * bw
  densN <- max(512, 2 * 2^(ceiling(log2(diff(rtRange)/(bw/2)))))
  endIdx <- 0
  maxFeatures <-  50
  num <- 0
  ###parallel processing of the density estination to speed up the process.
  vapp <- bpmapply((masspos[1:(length(masspos)-2)])[sel_idx],
                   (masspos[3:length(masspos)])[sel_idx]-1,
                   FUN=function(startIdx,endIdx,peaktable,
                                bw,maxFeatures,densFrom,densTo,densN){
                     if ((endIdx - startIdx) < 0)
                       return(matrix(0,nrow=0,ncol=2))

                     curMat <- peaktable[startIdx:endIdx, , drop = FALSE]

                     ###We generate the density for the model
                     den <- density(curMat[, 2], bw = bw, from = densFrom,
                                    to = densTo, n = densN,weights = curMat[,3])

                     deny <- den$y
                     maxy <- NULL
                     maxden <- max(deny)
                     snum <- 0
                     sub_res <- matrix(0,nrow=0,ncol=2)
                     vres <- vector(mode="list",length=maxFeatures)
                     for(i in seq_along(vres)){
                       vres[[i]] <- numeric(0)
                     }
                     num_feat <- 1
                     while (deny[maxy <- which.max(deny)] > maxden/20 && snum <
                            maxFeatures) {
                       if(num_feat > maxFeatures) break
                       grange <- xcms:::descendMin(deny, maxy)
                       deny[grange[1]:grange[2]] <- 0
                       peaktable_idx <- which(curMat[, 2] >= den$x[grange[1]] &
                                                curMat[, 2] <= den$x[grange[2]])
                       snum <- snum + 1
                       num <- num + 1
                       vres[[num_feat]] <- c(startIdx+peaktable_idx-1)
                       num_feat <- num_feat+1
                     }
                     if(num_feat>1){
                       return(vres[1:(num_feat-1)])
                     }
                    return(list())
                   },MoreArgs=list(peaktable=peaktable,
                                   bw=bw,maxFeatures=maxFeatures,
                                   densFrom=densFrom,densTo=densTo,densN),
                   BPPARAM = bpp,SIMPLIFY = FALSE)
  vapp <- do.call(c,vapp)
  vsize <- sapply(vapp,length)

  ###We retransform thge clustering table
  vapp <- mapply(vapp,seq_along(vapp),FUN = function(x,y){
    cbind(x,rep(y,length(x)))
  },SIMPLIFY = FALSE)

  vapp <- do.call(rbind,vapp)

  ###We select the best performing value
  vapp <- bplapply(split(as.data.frame(vapp),f = vapp[,1],
                 drop = FALSE),FUN=function(x,vsize){x[which.max(vsize[x[,2]]),]},
           vsize=vsize,BPPARAM = bpp)

  vapp <- do.call(rbind,vapp)

  groupindex <- rep(NA,nrow(peaktable))
  ###We convert the peaks back to the standard vector.
  vapp[,1] <- lam@order_peaks[vapp[,1]]
  vmm <- match(vapp[,2],unique(vapp[,2]))
  groupindex[vapp[,1]] <- vmm

  ###We corrct the missing value sif there is some
  vna <- which(is.na(groupindex))
  if(length(vna)>0){
    groupindex[vna] <- max(groupindex,na.rm = TRUE)+seq_along(vna)
  }
  return(groupindex)
  # lam@clustering <- groupindex
  ###We just return the index
  # return(lam)
}


###We correct the pealtable using the final clustering.
clusterPeaktable <- function(lam,clustering,bpp=NULL){
  temp_peaks <- by(lam@peaks,INDICES = clustering,FUN=function(x){
    c(sum(x[,1]*x[,3])/sum(x[,3]),sum(x[,2]*x[,3])/sum(x[,3]),sum(x[,3]),sum(x[,4]))
  })

  # new_idx <- match(1:nrow(lam@peaks),lam@order_peaks)

  temp_peaks <- do.call(rbind,temp_peaks)
  temp_peaks <- cbind(temp_peaks,1:nrow(temp_peaks))
  colnames(temp_peaks) <- c("mz","rt","num","current_num","id")
  lam@peaks <-as.data.frame(temp_peaks)

  if(is.null(bpp)) bpp <- bpparam()
  ###We correct all the blocks files.
  block_files <- getBlocksFiles(lam@storage, 1:length(lam@files))

  bplapply(as.list(block_files),FUN = function(x,index){
    tab <- read.table(x,header=TRUE,sep = ",")
    tab$id <- index[tab$id]
    write.table(tab,
                file = x,
                row.names = FALSE,
                sep = ",")
    return(NULL)
  },index=clustering,BPPARAM = bpp)
  return(lam)
}


####We eventually add the peaktable to the model.
addSinglePeaksToModel <- function(lam, peaks, cor_peaks, id_sample) {
  ###Extension of data
  npeaks <-
    data.frame(
      mz = peaks[, 1],
      rt = peaks[, 2],
      rtdev = rep(0, nrow(peaks)),
      rtfactor = rep(1, nrow(peaks)),
      num = rep(1, nrow(peaks)),current_num = rep(1,nrow(peaks)),
      id= 1:nrow(peaks)+lam@max_id,
      stringsAsFactors = FALSE)

  lam@max_id <- lam@max_id+nrow(peaks)
  if(nrow(lam@peaks)!=0){
    lam@peaks <- rbind(lam@peaks, npeaks)
  }else{
    lam@peaks <- npeaks
  }

  ###WE directly add the new feature to the table
  supp_data <- vector(mode="list",length=nrow(peaks))
  for(i in 1:length(supp_data)){
    x <- peaks[i,]
    temp <- data.frame(mz=c(x[[1]],0.0),rt=c(x[[2]],0.0),
                       rt_cor=c(cor_peaks[i,2],0.0), int=c(x[[3]],0.0),
                       sample=c(id_sample,NA_real_),stringsAsFactors = FALSE)
    supp_data[[i]] <-  temp
  }
  lam@data <- c(lam@data, supp_data)
  return(lam)
}

addPeaktableToModel <- function(lam, peaktable, cor_peaktable, index, id_sample, supp_data=character(0)) {
  ###We add the solo hypothesis on the data
  psolo <- is.na(index)
  to_insert <- which(!psolo)
  message("Found ",length(to_insert)," peaks belonging to the model and ",sum(psolo)," new peaks")

  ###Matche pdeak are added to the model corrected and uncorrected.
  if (length(to_insert) != 0){

    lam <- updateFeatureInfos(lam, cor_peaktable[to_insert, , drop = FALSE], index[to_insert])
  }

  ####New feature are directly added to the 2 matrices
  if(sum(psolo)!=0){
    psolo <- which(psolo)
    new_ids <- 1:length(psolo)+lam@max_id
    npeaks <-
      data.frame(
        mz = peaktable[psolo, 1],
        rt = cor_peaktable[psolo, 2],
        num = rep(1, length(psolo)),current_num = rep(1,length(psolo)),
        id= new_ids,
        stringsAsFactors = FALSE)
    lam@max_id <- lam@max_id+length(psolo)
    if(nrow(lam@peaks)!=0){
      lam@peaks <- rbind(lam@peaks, npeaks)
    }else{
      lam@peaks <- npeaks
    }
    index[psolo] <- new_ids
  }


  ###We habdle the case of supplementary informations evneutally.
  ldata <- list(peaktable[,1], peaktable[,2],
  cor_peaktable[,2],peaktable[,3])

  for(i in seq_along(supp_data)){
    ldata[[length(ldata)+1]] <- peaktable[,supp_data[i]]
  }
  ldata[[length(ldata)+1]] <- rep(id_sample,length(index))
  ldata[[length(ldata)+1]] <- index
  ndata <- do.call(cbind,ldata)

  colnames(ndata) <-   c('mz','rt','rt_cor','int',supp_data,'sample','id')
  lam@data <- rbind(lam@data,ndata)

  return(lam)
}




createFeatureInfos <- function(peaks, pos) {
  vp <- by(cbind(peaks,pos),INDICES = pos,FUN = function(x){
    return(c(mean(x[,1]),mean(x[,2]),nrow(x),nrow(x),x[1,ncol(x)]))
  })
  vp <- do.call(rbind,vp)
  colnames(vp)<- c("mz","rt",
                   "num","current_num","id")
  vp <- as.data.frame(vp)
  return(vp)
}



addPeaktablesToModel<- function(lam, peaktables, cor_peaktables, index, samples ,supp_data=character(0)){
  ###SIngle peaks
  psolo <- is.na(index)
  to_insert <- which(!psolo)
  
  ##New sample only features
  to_insert <- which((!is.na(index))&(index>nrow(lam@peaks)))
  
  ###Matched features
  to_update <- which((!is.na(index))&(index<=nrow(lam@peaks)))
  
  message("Found ",length(to_update)," peaks belonging to the model and ",sum(psolo)+length(unique(to_insert))," new features")
  
  ###Matche pdeak are added to the model corrected and uncorrected.
  if (length(to_update) > 0){
    
    lam <- updateFeatureInfos(lam, cor_peaktables[to_update, , drop = FALSE], index[to_update])
  }
  
  ###We aggregate the information of the new sample only peaks.
  if(length(to_insert)>0){
    ninfos <- createFeatureInfos(cor_peaktables[to_insert, , drop = FALSE], index[to_insert]) 
    lam@peaks <- rbind(lam@peaks, ninfos)
    lam@max_id <- lam@max_id+nrow(ninfos)
  }
  
  ####New feature are directly added to the 2 matrices
  if(sum(psolo)!=0){
    psolo <- which(psolo)
    new_ids <- 1:length(psolo)+lam@max_id
    npeaks <-
      data.frame(
        mz = peaktables[psolo, 1],
        rt = cor_peaktables[psolo, 2],
        num = rep(1, length(psolo)),current_num = rep(1,length(psolo)),
        id= new_ids,
        stringsAsFactors = FALSE)
    lam@max_id <- lam@max_id+length(psolo)
    if(nrow(lam@peaks)!=0){
      lam@peaks <- rbind(lam@peaks, npeaks)
    }else{
      lam@peaks <- npeaks
    }
    index[psolo] <- new_ids
  }
  

  
  ###We habdle the case of supplementary informations evneutally.
  ldata <- list(peaktables[,1], peaktables[,2],
                cor_peaktables[,2],peaktables[,3])
  
  for(i in seq_along(supp_data)){
    ldata[[length(ldata)+1]] <- peaktables[,supp_data[i]]
  }
  ldata[[length(ldata)+1]] <- samples
  ldata[[length(ldata)+1]] <- index
  ndata <- do.call(cbind,ldata)
  
  colnames(ndata) <-   c('mz','rt','rt_cor','int',supp_data,'sample','id')
  lam@data <- rbind(lam@data,ndata)
  
  return(lam)
}

setMethod("addFile","LCMSAlignerModel",function(object,path) {
  temp <- addFile(object@files, path = path)
  added <- FALSE
  if (temp$added) {
    object@files <- temp$obj
    added <- TRUE
  }
  return(list(obj = object, added = added, id = temp$id))
})

###WE need to store the RTs and the corrected RTs
#' Title
#'
#' @param lam An LCMS aligner object
#' @param path The path of the peaktableot be aligned
#' @param alignment_threshold The threshold used to align peaks
#' @param rt_scaling The RT scaling factor
#' @param ransac_niter The number of iteration of the RANSAC
#' @param ransac_dist_threshold The distance trheshold used for RNASAC inlier detection
#' @param rt_extensions The retention time of RT
#'
#' @return
#' @export
#'
#' @examples
#' print('Examples to be put here')
alignPeaktable <-
  function(lam,
           path,
           rt_scaling = c(0.9, 0.95, 1, 1.05, 1.1),
           ransac_niter = 2000,
           ransac_dist_threshold = NULL,
           ransac_l1=0.05,
           rt_extensions = c(4, 0.15),
           path_aligner=NULL,span = 0.6,
           supp_data=c("peakwidth","SN","right_on_left_assymetry"),
           lim = 0.1,
           graphical=FALSE,
           bpp=NULL) {
    ##Case of first added file.
    temp <- addFile(lam, path)
    if (temp$added) {
      lam <- temp$obj
    }else{
      return(lam)
    }

    if(is.null(ransac_dist_threshold)) ransac_dist_threshold <- lam@references@parameters$rt/2
    lam@files_in_memory <- c(lam@files_in_memory,temp$id)
    alignment_threshold <- lam@references@parameters[["alignment_threshold"]]

    ###WE retrieve the id associated to the sample
    id <- temp$id
    message("Aligning ",path)
    peaktable <- read.table(path,sep = lam@references@parameters$sep_table,header=TRUE)
    peaktable <- fortifyPeaktable(peaktable,lam@references,supp_data=supp_data)
    ###We first correct the retention time
     lcc <-
      correctPeaktable(
        lam@references,
        peaktable,
        lam@peaks,
        rt_scaling = rt_scaling,
        ransac_niter = ransac_niter,
        ransac_dist_threshold = ransac_dist_threshold,
        ransac_l1=ransac_l1,
        extensions = rt_extensions,
        lim = lim
      )
    cpeaktable <- lcc$peaktable
    ###We reorder the two peaktables
    opeaktable <- order(cpeaktable[,"mz"])
    cpeaktable <- cpeaktable[opeaktable,,drop=FALSE]
    peaktable <- peaktable[opeaktable,,drop=FALSE]

    aoo <- order(peaktable[,2])
    vsel <- ceiling(seq(1,length(aoo)-1,length=150))
    tcor <- matrix(c(peaktable[aoo[vsel],2],cpeaktable[aoo[vsel],2]-peaktable[aoo[vsel],2]),ncol=2,byrow=FALSE)
    lam@rt_correction[[length(lam@files)]] <- tcor

    vmatch <- NULL

    ###We always reorder the peaks by mass

    if(nrow(lam@peaks)!=0){
      if(lam@references@parameters$aligner_type=="density"){
        #####The bwe is always 2 times the retention time ost
        vmatch <- alignToModel.density(lam, cpeaktable,bw = lam@references@parameters$rt_dens,
                                       binSize = lam@references@parameters$dmz,
                                       maxFeatures = 50, sleep = 0,
                                       mzCost = lam@references@parameters$dmz,
                                       rtCost = lam@references@parameters$rt, bpp = bpp)
      }else if(lam@references@parameters$aligner_type=="nn"){
        vmatch <-
          alignToModel.nn(lam, cpeaktable, threshold = alignment_threshold,id_sample=id)
      }
    }else{
      vmatch <- rep(NA_real_,nrow(peaktable))
    }
    ###We can then change the match vector to the non sorted coordinates
    pfound <- !is.na(vmatch)
    vmatch[pfound] <- lam@order_peaks[vmatch[pfound]]

    ##We check if there is a major difference between the mapped feature and the rest
    if(any(pfound)){
      tdm <- cbind(lam@peaks[vmatch[pfound],1],cpeaktable[pfound,1])
      tdm <- apply(tdm,1,function(x){diff(range(x))})
    }
    
    lam <- addPeaktableToModel(lam, peaktable, cpeaktable, vmatch,id_sample = id,supp_data=supp_data)

    lam@order_peaks <- order(lam@peaks[,1])
    if(!is.null(path_aligner)&&(length(lam@files_in_memory)>=lam@save_interval)){
      lam <- saveAligner(lam,path_aligner,supp_data=supp_data,reset=TRUE)
    }
    # message("Peaktable row:",nrow(cpeaktable)," peaks model:",nrow(lam@peaks),
    #         "total num",sum(lam@peaks[,"num"]))

    return(lam)
  }



#' Title
#'
#' @param lam An LCMS aligner object
#' @param path The path of the peaktableot be aligned
#' @param alignment_threshold The threshold used to align peaks
#' @param rt_scaling The RT scaling factor
#' @param ransac_niter The number of iteration of the RANSAC
#' @param ransac_dist_threshold The distance trheshold used for RNASAC inlier detection
#' @param rt_extensions The retention time of RT
#'
#' @return
#' @export
#'
#' @examples
#' print('Examples to be put here')
alignPeaktables <-
  function(lam,
           paths,
           rt_scaling = c(0.9, 0.95, 1, 1.05, 1.1),
           ransac_niter = 1500,
           ransac_dist_threshold = NULL,
           ransac_l1=0.05,
           ransac_span=0.7,
           rt_extensions = c(4, 0.15),
           path_aligner=NULL,span = 0.6,
           supp_data=c("peakwidth","SN","right_on_left_assymetry"),
           lim = 0.1,
           graphical=FALSE,
           correct_rt=TRUE,
           bpp=NULL) {
    message("Starting batch.")
    start_time <- Sys.time()
    paths <- paths[!is.na(paths)]
    if(length(paths)==0) return(lam)
    ##Case of first added file.
    added_paths <- rep(TRUE,length(paths))
    ids <- rep(NA,length(paths))
    for(ipath in seq_along(paths)){
      path <- paths[ipath]
      temp <- addFile(lam, path)
      if (temp$added) {
        ids[ipath] <- temp$id
        lam <- temp$obj
      }else{
        added_paths[ipath] <- FALSE
      }
    }
    
    if(all(!added_paths)) return(lam)
    ###We remove the previousoly added files
    paths <- paths[added_paths]
    ids <- ids[added_paths]
    
    if(is.null(ransac_dist_threshold)) ransac_dist_threshold <- lam@references@parameters$rt/2
    lam@files_in_memory <- c(lam@files_in_memory,ids)
    alignment_threshold <- lam@references@parameters[["alignment_threshold"]]
    
    # message("Aligning ",length(paths)," files. CUrrent aligner size is ",format(object.size(lam),"Mb"))
    
    correctPeaktablePar <- function(path,ref,supp_data,rt_scaling,ransac_niter,
                                 ransac_dist_threshold,ransac_l1,rt_extensions,ransac_span,
                                 lim,graphical,peaks,ratio,correct_rt){
      if(!exists("fortifyPeaktable")){
        library(onlineLCMSaligner)
      }
      peaktable <- read.table(path,sep = ref@parameters$sep_table,header=TRUE)
      if(nrow(peaktable)==0){
        return(list(peaktable,peaktable,matrix(c(NA,NA),ncol=2,nrow=1)))
        
      }
      
      # message("Correctign file: ",basename(path)," against ",nrow(peaks)," features.")
      
      peaktable <- fortifyPeaktable(peaktable,ref,supp_data=supp_data)
      
      cpeaktable <- NULL
      if(correct_rt & nrow(peaktable)>2){
      ###We first correct the retention time
      lcc <-
        tryCatch(correctPeaktable(
          ref,
          peaktable,
          peaks=peaks,
          rt_scaling = rt_scaling,
          ransac_niter = ransac_niter,
          ransac_dist_threshold = ransac_dist_threshold,
          ransac_l1=ransac_l1,
          ransac_span=ransac_span,
          extensions = rt_extensions,
          ratio = ratio,
          lim = lim,
          graphical = graphical
        ),error=function(e){return(NA)})
        if(is.na(lcc)){
          cpeaktable <- peaktable
        }else{
          cpeaktable <- lcc$peaktable
        }
      }else{
        cpeaktable <- peaktable
      }
      
      ###We reorder the two peaktables
      opeaktable <- order(cpeaktable[,"mz"])
      cpeaktable <- cpeaktable[opeaktable,,drop=FALSE]
      peaktable <- peaktable[opeaktable,,drop=FALSE]
      
      aoo <- order(peaktable[,2])
      if(length(aoo)<=10){
        tcor <- matrix(c(peaktable[aoo,2],cpeaktable[,2]-peaktable[aoo,2]),ncol=2,byrow=FALSE)
        return(list(peaktable,cpeaktable,tcor))
      }
      vsel <- ceiling(seq(1,length(aoo)-1,length=150))
      tcor <- matrix(c(peaktable[aoo[vsel],2],cpeaktable[aoo[vsel],2]-peaktable[aoo[vsel],2]),ncol=2,byrow=FALSE)
      return(list(peaktable,cpeaktable,tcor))
    }
    tpeaks <- lam@peaks
    ratio <- 1
    if(nrow(tpeaks)>10000){
      ratio <- 10000/nrow(tpeaks)
      tpeaks <- tpeaks[order(tpeaks$num,decreasing=TRUE)[1:5000],]
    }
    
    message("Starting retention correction.")
    start_time <- Sys.time()
    values <- bplapply(paths,FUN = correctPeaktablePar,ref=lam@references,rt_scaling = rt_scaling,
             ransac_niter = ransac_niter,
             ransac_dist_threshold = ransac_dist_threshold,
             ransac_l1=ransac_l1,peaks=tpeaks,
             rt_extensions = rt_extensions, ratio=ratio,
             ransac_span = ransac_span,lim = lim,
             supp_data=supp_data,graphical=graphical,correct_rt=correct_rt,BPPARAM = bpp)
    end_time <- Sys.time()
    # message("Retention time correction finshed in ", end_time-start_time))
    ### We bind all the peaktable
    peaktables <- lapply(values,"[[",1)
    
    ###If it is the first peak we use it as a reference.
    if(nrow(lam@peaks)==1){
      idref <- 1
      while(nrow(peaktables[[idref]])==0){
        idref <- idref+1
        ##All peaktables are empty
        if(idref>length(peaktables)) return(lam)
      }
      ##We extract the first non empty peaktable
      pt_ref <- peaktables[[idref]]
      peaktables <- peaktables[-idref]
      vpt_ref <- cpeaktables[[idref]]
      cpeaktables <- cpeaktables[-idref]
      
      lam <- addPeaktableToModel(lam, pt_ref, cpt_ref, rep(NA,1:nrow(pt_ref)), 1, supp_data=supp_data)
    }
    
    
    
    sample_bins <- sapply(peaktables,nrow)
    sample_vec <- rep(ids,times=sample_bins)
    
    
    peaktables <- do.call(rbind,peaktables)
    
    cpeaktables <-  lapply(values,"[[",2)
    cpeaktables <- do.call(rbind,cpeaktables)
    rt_corr <- lapply(values,"[[",3)
    
    
    lam@rt_correction[(length(lam@files)-length(paths)+1):length(lam@files)] <- rt_corr
    
    
    ###We reorder the two peaktables
    opeaktables <- order(cpeaktables[,"mz"])
    cpeaktables <- cpeaktables[opeaktables,,drop=FALSE]
    peaktables <- peaktables[opeaktables,,drop=FALSE]
    sample_vec <- sample_vec[opeaktables]


    vmatch <- NULL
    ###We match the features
    
    #####The bwe is always 2 times the retention time ost
    message("Computing density")
    start_time <- Sys.time()
    vmatch <- alignToModelByBatch.density(lam, cpeaktables,bw = lam@references@parameters$rt_dens,
                                   binSize = lam@references@parameters$dmz,
                                   maxFeatures = 50, sleep = 0,
                                   mzCost = lam@references@parameters$dmz,
                                   rtCost = lam@references@parameters$rt, bpp = bpp)
    
    end_time <- Sys.time()
    ###We can then change the match vector to the non sorted order
    pfound <- which(!is.na(vmatch))
    pfound <- pfound[vmatch[pfound]<=nrow(lam@peaks)]
    if(any(pfound)){
      vmatch[pfound] <- lam@order_peaks[vmatch[pfound]]
    }
    message("Adding aligned features to the model.")
    message("Starting adding aligned features to the model ",Sys.time())
    lam <- addPeaktablesToModel(lam, peaktables, cpeaktables, vmatch,
                                sample_vec,supp_data=supp_data)
    message("Ending adding aligned features to the model ",Sys.time())
    lam@order_peaks <- order(lam@peaks[,1])
    if(!is.null(path_aligner)&&(length(lam@files_in_memory)>=lam@save_interval)){
      lam <- saveAligner(lam,path_aligner,supp_data=supp_data,reset=TRUE)
    }
    # message("Peaktable row:",nrow(cpeaktable)," peaks model:",nrow(lam@peaks),
    #         "total num",sum(lam@peaks[,"num"]))
    
    return(lam)
  }


#' Writing an alignment to a file
#'
#' @param lam An LCMSAlignerModel
#' @param path The path ofhte alignment object to read
#'
#' @return A boolean indicating if the file has been succesfully written
#' @export
#'
#' @examples
#' print("Examples ot be put here")
writeAlignment <- function(lam, path) {
  saveRDS(lam, path)
  return(file.exists(path))
}

#' Reading an alignment from a file
#'
#' @param path The path ofhte alignment object to read
#'
#' @return An LCMSAlignerModel object
#' @export
#'
#' @examples
#' print("Examples ot be put here")
readAlignment <- function(path) {
  readRDS(path)
}


saveAligner <- function(lam,path,reset=TRUE,supp_data=character(0)){
  message("Saving alignment.")
  ####WE get the data to_write
  ####Write the aligned data.
  if(nrow(lam@data)>0){
    lam@storage <- writeBlock(lam@storage,lam@data,lam@files_in_memory)
  }
  ###We write the state of the aligner after writing the files
  ###Reset the data
  if(reset){
    lam <- reset(lam,supp_data=supp_data)
  }
  writeAlignment(lam,path)
  return(lam)
}

###Function to clean after the data have been written or
###because the processing needs to be restarted.
restart <- function(lam,supp_data=character(0)){
  ###We erase the storage eventually
  resetStorage(lam@storage)
  lam@storage <- LCMSAlignerDataStorage(lam@storage@dir)
  lam@save_interval <- lam@save_interval
  lam@peaks <- data.frame(mz = numeric(0), rt = numeric(0), 
                          num=numeric(0),current_num=numeric(0),
                          id = character(0),stringsAsFactors = FALSE)
  temp <- matrix(0,nrow=0,ncol=6)
  lam@data <- matrix(ncol = 6+length(supp_data), nrow = 0)
  colnames(lam@data) <- c('mz', 'rt', 'rt_cor', 'int',supp_data, 'sample', 'id')
  lam@data <- temp
  ###Input processing
  lam@files <- LCMSAlignerFileHandler()
  return(lam)
}


reset <- function(lam,supp_data=character(0)){
  ###We empty all the data matrix.
  temp <- matrix(0,nrow=0,ncol=6+length(supp_data))
  colnames(temp) <- c('mz','rt','rt_cor','int',supp_data,'sample','id')
  lam@data <- temp
  ###We reset the counter of peaks.
  for(i in seq_along(lam@data)){
    lam@peaks[i,"current_num"] <- 0
  }
  lam@files_in_memory <- numeric(0)
  return(lam)
}

setMethod("isProcessed","LCMSAlignerModel",function(object,paths){
  return(isProcessed(object@files,paths))
})



make_raw_name <- function(x){
  str_split(basename(x),pattern = fixed("."))[[1]][1]
}

####Data matrix building function

#' Title
#'
#' @param LCMSAlignerModel
#'
#' @return
#' @export
#'
#' @examples
#' print("example to be put here")
exportDataMatrix <- function(object,path,subsample=NULL,subvariable=NULL,
                                                        quant_var="int", add_summary=TRUE,
                                                        summary_vars= c("peakwidth","rt","rt_cor",
                                                                        "right_on_left_assymetry"),max_memory=2000){
  if(is.null(subsample)){
    subsample <- 1:length(object@files)
  }
  fnames <- getAllFiles(object@files)
  if(!is.numeric(subsample)){
    ###Give the position of the match on the data matrix
    subsample <- match(basename(subsample),basename(fnames))
  }
  mzs <- NULL
  rts <-  NULL
  if(is.null(subvariable)){
    subvariable <- 1:nrow(object@peaks)
  }
  if(is.logical(subvariable)){
    subvariable <- which(subvariable)
  }
  
    #pf <- match(subvariable,object@peaks$id)
  mzs <- object@peaks$mz[subvariable]
  rts <- object@peaks$rt[subvariable]
  
  all_order <- order(mzs,rts,decreasing = FALSE)
  
  names_cols <- paste(quant_var,basename(fnames),sep="_")
  names_cols <- names_cols[subsample]
  
  ###We split the data matrix to avoid filling the memoery too fast
  ##matrix size = 216+8*num_elemnts bytes
  BATCH_SIZE <- 20000
  chunks <- seq(1,length(mzs),by=BATCH_SIZE)
  if(chunks[length(chunks)]!=length(mzs)){
    chunks <- c(chunks,length(mzs))
  }
  chunks[length(chunks)] <- chunks[length(chunks)]+1
  
  for(ic in 1:(length(chunks)-1)){
    sub_idx <- chunks[ic]:(chunks[ic+1]-1)
  
    dmat <- buildDataMatrix.datatable(object@storage,subsample=subsample,
                            subvariable=subvariable[all_order][sub_idx],max_sample=length(object@files),
                            quant_var=quant_var,summary_vars=summary_vars,name_samples = names_cols,
                            bpp =object@references@bpp)
    # colnames(dmat) <- names_cols
    cnames <- colnames(dmat)
    dmat <- cbind(mzs[all_order][sub_idx],rts[all_order][sub_idx],dmat)
    colnames(dmat) <- c("mz","rt",cnames)
    ###We know write the data
    if(ic==1){
      write_dgCMatrix_csv(dmat, path,chunk_size = 2000,replace_old = 0.0,replace_new = NA, append = FALSE)
    }else{
      write_dgCMatrix_csv(dmat, path,chunk_size = 2000,replace_old = 0.0,replace_new = NA, append = TRUE)
    }
  }
}



#' The percentage of alignement in the data
#'
#' @param lam The LCMSalignerModel to align
#' @param int_threshold The alignement threshold to be used eventually.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' print("Examples to be put here")
plotDevRt <- function(lam,int_threshold= 0.15){
  maxr <- sapply(lam@rt_correction,function(x){return(c(range(x[,1]),range(x[,2])))})
  mint <- min(maxr[1,],na.rm = TRUE)
  maxt <- max(maxr[2,],na.rm = TRUE)
  mindt <- min(maxr[3,],na.rm = TRUE)
  maxdt <- max(maxr[4,],na.rm = TRUE)
  meandt <- maxdt-mindt

  ylim <- NULL
  mdev <- int_threshold*(maxt-mint)
  if(max(meandt-mindt,maxdt-meandt)<mdev){
    ylim <- c(meandt-mdev,meandt+mdev)
  }else{
    ylim <- c(mindt,maxdt)
  }

  plot(NULL,xlim=c(mint,maxt),ylim=ylim,xlab="RT(mins)",ylab="dRT(mins)",main="Retention time correction")
  colv <- rainbow(length(lam@rt_correction))
  for(i in seq_along(colv)){
    if(is.na(lam@rt_correction[[i]][1,1])) next
    lines(lam@rt_correction[[i]][,1],lam@rt_correction[[i]][,2],col=colv[i])
  }
}
