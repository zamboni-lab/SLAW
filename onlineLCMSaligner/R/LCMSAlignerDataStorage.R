#' @include classes.R
#' @include allGenerics.R
#' @include utils.R

LCMSAlignerDataStorage <- function(path, force = TRUE) {
  if (!dir.exists(path)) {
    warning("path directory doe snot exists, it is created.")
    dir.create(path)
  } else{
    allf <- list.files()

  }
  temp <-  new("LCMSAlignerDataStorage")
  temp@dir <- path
  temp@files <- list()
  return(temp)
}

# setClass("LCMSAlignerDataStorage",
#          slots = list(dir="character",files="list",frequency="numeric"))

storeData <- function(data, filename) {
  ####We write the data to a file
  fwrite(as.data.table(data),
              file = filename,
              row.names = FALSE,
              sep = ",")
  # message("data stored")
}

makeFileName <- function(lads) {
  fname <-
    file.path(lads@dir, paste("block_", length(lads@files) + 1, ".csv", sep =
                                ""))
  return(fname)
}


setMethod("length", "LCMSAlignerDataStorage", function(x) {
  length(x@files)

})

####Write a new block of peaks to the disk
writeBlock <- function(lads, data, id_files) {
  new_name <- makeFileName(lads)
  # message("name made")
  lads@files[[new_name]] <- id_files
  sdf <- storeData(data, new_name)
  # message("block written")
  ###Normally this is enough to sample all the data
  return(lads)
}


###In this section subsample is always numeric

extractSubDataMatrix.data.table <- function(path, subsample, subvariable) {
  ###We only keep the
  vtable <-
    fread(path,
               header = TRUE,
               sep = ",")
  
  psel <-
    (vtable[, "id"] %in% subvariable) &
    (vtable[, "sample"] %in% subsample)
  
  return(vtable[sample %in% subsample & id %in% subvariable,])
}

extractSubDataMatrix <- function(path, subsample, subvariable) {
  ###We only keep the
  vtable <-
    read.table(path,
               header = TRUE,
               sep = ",",
               stringsAsFactors = FALSE)
  psel <-
    (vtable[, "id"] %in% subvariable) &
    (vtable[, "sample"] %in% subsample)
  return(vtable[psel, , drop = FALSE])
}


getBlocksFiles <- function(lads, subsample) {
  sel_blocks <- sapply(lads@files, function(x, ref) {
    any(x %in% ref)
  }, ref = subsample)
  return(names(lads@files)[sel_blocks])
}


resetStorage <- function(lads) {
  lfa <- list.files(lads@dir, full.names = TRUE)
  if (length(lfa) > 0)
    do.call(file.remove, as.list(lfa))
  return(lads)
}

###Subsample is jsut a set of integer corresponding to the idx of the smaple let subsample tpo 0 eventually
###subvariable corresponds to the selected variable in the right order
 ###We also summarize the submatrix eventually.
setMethod("buildDataMatrix", "LCMSAlignerDataStorage", function(object,
                                                                subsample,
                                                                subvariable,
                                                                summary_vars,
                                                                quant_var,
                                                                name_samples,
                                                                bpp=NULL,
                                                                max_sample = max(subsample)
                                                                ){
  if(is.null(bpp)) bpp <- bpparam()
  ###number and name of the file are in the same order
  ###The subsample file is always extracted of the data
  new_idx_sample <- rep(NA_real_, length(subsample))
  new_idx_sample[subsample] <- seq_along(subsample)

  ###We define the new indices of the features matrix
  new_idx_variables <- sort(subvariable, decreasing = FALSE)

  ###We get a list of the blocks file
  block_files <- getBlocksFiles(object, subsample)


  ###We remove the rest of the values
  block <- extractSubDataMatrix(block_files[1], subsample, subvariable)
  sel_vars <- match(summary_vars,colnames(block))
  sel_pos <- which(!is.na(sel_vars))
  sel_order <- seq_along(sel_pos)
  sel_blocks <- sel_vars[sel_pos]

  ###THis is to store all the data in this.
  mu <- matrix(0,nrow=length(subvariable),ncol=length(sel_blocks))
  vmin <- matrix(1000000,nrow=length(subvariable),ncol=length(sel_blocks))
  vmax <- matrix(0,nrow=length(subvariable),ncol=length(sel_blocks))
  nums <- rep(0,length=length(subvariable))
  # sds <- matrix(0,nrow=length(subvariable),ncol=length(sel_blocks))


# 
#   if ((memory.limit() * 10 ^ 6 * 0.8) < (length(subvariable) * (length(subsample)+3*length(summary_vars)) *
#                                          8)) {
#     stop("Matrix to be written probably too big to fit in memoery, please filter some variable.")
#   }

  ###Columns are mz rt and rt_dev among acquisition
  vmat <-
    matrix(NA_real_,
           nrow = length(subvariable),
           ncol = length(subsample))


  id_mean <- 3:(length(sel_blocks)+2)
  id_min <- (length(sel_blocks)+3):(2*length(sel_blocks)+2)
  id_max <- (2*length(sel_blocks)+3):(3*length(sel_blocks)+2)

  ###We add m/z and the retention time
  for (ff in block_files) {
    ###We read the block
    block <- extractSubDataMatrix(ff, subsample, subvariable)

    ###updateing the intensites
    vagg_int <- bplapply(split(block,f =paste(block[, "id"],block[,"sample"],sep="_"),drop = FALSE),
                         FUN = function(x,var_quant){
      c(max(x[,var_quant]),x[1,"id"],x[1,"sample"])
    },var_quant=quant_var,BPPARAM = bpp)

    # vagg_int <- by(data = block,FUN = function(x,var_quant){
    #   c(max(x[,var_quant]),x[1,"id"],x[1,"sample"])
    # },INDICES=paste(block[, "id"],block[,"sample"],sep="_"),var_quant=quant_var)
    vagg_int <- do.call(rbind,vagg_int)

    idcol <- match(vagg_int[,3], new_idx_sample)
    idrow <- match(vagg_int[,2], new_idx_variables)

    ###We update al the intensities
    idx_mat <- matrix(c(idrow,idcol),ncol=2,byrow=FALSE)
    vmat[idx_mat] <- apply(cbind(vmat[idx_mat],vagg_int[,1]),1,max,na.rm=TRUE)

    ###updating all the other variables
    ###we create a summary table for faster computation.
    # vagg <- by(data = block,FUN = function(x,subvars){
    #   c(length(x),x[1,"id"],apply(x[,subvars],2,mean),apply(x[,subvars],2,min),
    #     apply(x[,subvars],2,max))
    # },INDICES=block[, "id"],subvars=sel_blocks)
    # browser()

    vagg <- bplapply(split(block,f = block[,"id"],drop = FALSE),FUN = function(x,subvars){
         c(nrow(x),x[1,"id"],apply(x[,subvars],2,mean),apply(x[,subvars],2,min),
              +       apply(x[,subvars],2,max))
      },subvars=sel_blocks,BPPARAM = bpp)

    vagg <- do.call(rbind,vagg)

    idrowb <- match(vagg[,2], new_idx_variables)

    ###We update all the other variables
    nums[idrowb] <- nums[idrowb]+vagg[,1]


    ###We update the mu estimate by summing only
    idx_mat_2 <- cbind(rep(idrowb,times=length(sel_order)),rep(sel_order,each=nrow(vagg)))
    mu[idx_mat_2] <- mu[idx_mat_2]+vagg[,1]*as.numeric(vagg[,id_mean])
    vmin[idx_mat_2] <- apply(cbind(vmin[idx_mat_2],as.numeric(vagg[,id_min])),1,min)
    vmax[idx_mat_2] <- apply(cbind(vmax[idx_mat_2],as.numeric(vagg[,id_max])),1,max)
  }

  # ###We update the mean calculations
  mu <- mu/nums
  #
  # ###second opening of file to dertermine variances.
  # for (ff in block_files) {
  #   block <- extractSubDataMatrix(ff, subsample, subvariable)
  #   idrow <- match(block[, "id"], new_idx_variables)
  #   for(idf in seq_along(idrow)){
  #     sds[idrow[idf],sel_order] <- as.numeric((block[idf,sel_blocks]-mu[idrow[idf],sel_order])**2+sds[idrow[idf],sel_order])
  #   }
  # }
  # sds <- sqrt(sds/nums)
  # cvs <- sds/mu
  colnames(vmat) <- name_samples
  rmat <- cbind(nums,mu,vmin,vmax)
  summary_vars <- summary_vars[sel_pos]
  cnames <- c("num_detection",paste("mean",summary_vars,sep="_"),
              paste("min",summary_vars,sep="_"),paste("max",summary_vars,sep="_"))
  ###We shuffle the integre sequence evneutally
  new_order <- c(1,as.numeric(t(matrix(1:(ncol(rmat)-1),ncol=3)))+1)
  colnames(rmat) <- cnames
  rmat <- rmat[,new_order,drop=FALSE]
  final_mat <- cbind(rmat,vmat)
  ###The number of detection is actualized.
  final_mat[,"num_detection"] <- apply(vmat,1,function(x){sum(!is.na(x))})
  return(final_mat)
})



###Subsample is jsut a set of integer corresponding to the idx of the smaple let subsample tpo 0 eventually
###subvariable corresponds to the selected variable in the right order
###We also summarize the submatrix eventually.




# TEST_PATH <- c("C:/Users/dalexis/Documents/dev/docker_sing_output/block_98.csv",
#                "C:/Users/dalexis/Documents/dev/docker_sing_output/block_133.csv")
# tt <- readRDS("C:/Users/dalexis/Documents/dev/docker_sing_output/alignement.rds")
# subvariable <- which(tt@peaks$num>1500)
# subsample <- 1:2400
# path <- TEST_PATH[1]

buildDataMatrix.datatable <- function(object,
                                                                subsample,
                                                                subvariable,
                                                                summary_vars,
                                                                quant_var,
                                                                name_samples,
                                                                bpp=NULL,
                                                                max_sample = max(subsample)
){
  if(is.null(bpp)) bpp <- bpparam()
  ###number and name of the file are in the same order
  ###The subsample file is always extracted of the data
  new_idx_sample <- rep(NA_real_, length(subsample))
  new_idx_sample[subsample] <- seq_along(subsample)
  
  ###We get a list of the blocks file
  block_files <- getBlocksFiles(object, subsample)
  
  
  ###We remove the rest of the values
  block <- extractSubDataMatrix.data.table(block_files[1], subsample, subvariable)
  sel_vars <- match(summary_vars,colnames(block))
  sel_pos <- which(!is.na(sel_vars))
  sel_order <- seq_along(sel_pos)
  sel_blocks <- sel_vars[sel_pos]
  
  ###THis is to store all the data in this.
  mu <- matrix(0,nrow=length(subvariable),ncol=length(sel_blocks))
  vmin <- matrix(1000000,nrow=length(subvariable),ncol=length(sel_blocks))
  vmax <- matrix(0,nrow=length(subvariable),ncol=length(sel_blocks))
  nums <- rep(0,length=length(subvariable))
  # sds <- matrix(0,nrow=length(subvariable),ncol=length(sel_blocks))
  
  
  # 
  # if ((memory.limit() * 10 ^ 6 * 0.8) < (length(subvariable) * (length(subsample)+3*length(summary_vars)) *
  #                                        8)) {
  #   stop("Matrix to be written probably too big to fit in memoery, please filter some variable.")
  # }
  # 
  ###Columns are mz rt and rt_dev among acquisition
  vmat <- Matrix(0.0,nrow = length(subvariable),
                 ncol = length(subsample),sparse = TRUE)
  # vm <- data.table(nrow=length(subvariable),ncol=length(subsample))
  
  
  id_mean <- 3:(length(sel_blocks)+2)
  id_min <- (length(sel_blocks)+3):(2*length(sel_blocks)+2)
  id_max <- (2*length(sel_blocks)+3):(3*length(sel_blocks)+2)
  
  ###We add m/z and the retention time
  for (ff in block_files) {
    ###We read the block
    block <- extractSubDataMatrix.data.table(ff, subsample, subvariable)
    tff <- function(x,var_quant){
      return(c(max(x[, ..var_quant])))
    }
    vagg_int <- block[, tff(.SD,quant_var), by = .(id, sample)]
    
    idcol <- match(vagg_int[,sample], new_idx_sample)
    idrow <- match(vagg_int[,id], subvariable)
    
    mat_idx <- cbind(idrow,idcol)
    ###We update al the intensities
    idx_mat <- matrix(c(idrow,idcol),ncol=2,byrow=FALSE)
    vmat[idx_mat] <- apply(cbind(vmat[idx_mat],vagg_int[[3]]),1,max,na.rm=TRUE)
    
    tff2 <- function(x,subvars){
      c(nrow(x),apply(x[,..subvars],2,mean),apply(x[,..subvars],2,min),apply(x[,..subvars],2,max))
    }
    
    vagg <- block[, as.list(tff2(.SD,sel_blocks)), by = .(id)]

    
    idrowb <- match(vagg[[1]], subvariable)
    
    ###We update all the other variables
    nums[idrowb] <- nums[idrowb]+vagg[[2]]
    
    
    ###We update the mu estimate by summing only
    idx_mat_2 <- cbind(rep(idrowb,times=length(sel_order)),rep(sel_order,each=nrow(vagg)))
    mu[idx_mat_2] <- mu[idx_mat_2]+vagg[[2]]*as.numeric(as.matrix(vagg[,..id_mean]))
    vmin[idx_mat_2] <- apply(cbind(vmin[idx_mat_2],as.numeric(as.matrix(vagg[,..id_min]))),1,min)
    vmax[idx_mat_2] <- apply(cbind(vmax[idx_mat_2],as.numeric(as.matrix(vagg[,..id_max]))),1,max)
  }
  
  # ###We update the mean calculations
  sel_mu <- mu[,1]!=0
  mu[sel_mu] <- mu[sel_mu]/nums[sel_mu]
  colnames(vmat) <- name_samples
  
  summary_vars <- summary_vars[sel_pos]
  unique_detection <- apply(vmat,1,function(x){sum(x!=0)})
  rmat <- cbind(unique_detection,nums,mu,vmin,vmax)
  cnames <- c("num_detection","total_detection",paste("mean",summary_vars,sep="_"),
              paste("min",summary_vars,sep="_"),paste("max",summary_vars,sep="_"))
  tsidx <- seq_along(summary_vars)
  new_order <- as.numeric(rbind(tsidx,tsidx+length(summary_vars),tsidx+2*length(summary_vars)))
  new_order <- c(1,2,new_order+2)
  colnames(rmat) <- cnames
  rmat <- rmat[,new_order,drop=FALSE]
  final_mat <- cbind(rmat,vmat)
  return(final_mat)
}




setMethod("summarizeMetrics","LCMSAlignerDataStorage",function(object,subsample,subvariable,vars){
  ###We calculate the mean and the standard eviation of all checked qunatitiy online.

  ###We define the new indices of the features matrix
  new_idx_variables <- sort(subvariable, decreasing = FALSE)

  ###We get a list of the blocks file
  block_files <- getBlocksFiles(object, subsample)

  ###We laways get the number of matched values
  block <- extractSubDataMatrix(block_files[1], subsample, subvariable)
  sel_vars <- match(vars,colnames(block))
  sel_pos <- which(!is.na(sel_vars))
  sel_order <- seq_along(sel_pos)
  sel_blocks <- sel_vars[sel_pos]

  ###Columns are mz rt and rt_dev among acquisition
  vmat <-
    matrix(NA_real_,
           nrow = length(subvariable),
           ncol = 3*length(sel_pos))



  ###We add m/z and the retention time
  for (ff in block_files) {
    ###We read the block
    block <- extractSubDataMatrix(ff, subsample, subvariable)
    idrow <- match(block[, "id"], new_idx_variables)

  }
  bool_pres <- FALSE
  ###second pass uting all the standard deviation.
  for (ff in block_files) {
    block <- extractSubDataMatrix(ff, subsample, subvariable)
    idrow <- match(block[, "id"], new_idx_variables)
    for(idf in seq_along(idrow)){
      sds[idrow[idf],sel_order] <- as.numeric((block[idf,sel_blocks]-mu[idrow[idf],sel_order])**2+sds[idrow[idf],sel_order])
    }
  }
  sds <- sqrt(sds/nums)
  cvs <- sds/mu

  rmat <- cbind(nums,mu,sds,cvs)
  vars <- vars[sel_pos]
  cnames <- c("num_detection",paste("mean",vars,sep="_"),paste("sd",vars,sep="_"),paste("cv",vars,sep="_"))
  ###We shuffle the integre sequence evneutally
  new_order <- c(1,as.numeric(t(matrix(1:(ncol(rmat)-1),ncol=3)))+1)
  colnames(rmat) <- cnames
  rmat <- rmat[,new_order,drop=FALSE]
})
