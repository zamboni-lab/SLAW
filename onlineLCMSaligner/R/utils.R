extendList <- function(ll,idx){
  if(length(ll)<idx){
    ll <- c(ll,vector(mode="list",length=length(ll)))
  }
  return(ll)
}

extendVector <- function(vv,idx){
  if(length(vv)<idx){
    vv <- c(vv,do.call(typeof(vv),list(length(vv))))
  }
  return(vv)
}


extendMatrixRow <- function(mm,idx){
  if(nrow(mm)<idx){
    mm <- rbind(mm,matrix(mm[1,1],ncol=ncol(mm),nrow=nrow(mm)))
  }
  return(mm)
}

extendDfRow <- function(mm,idx){
  if(nrow(mm)<idx){
    mm <- rbind(mm,mm)
  }
  return(mm)
}

extendMatrixCol <- function(mm,idx){
  if(ncol(mm)<idx){
    mm <- cbind(mm,matrix(mm[1,1],ncol=ncol(mm),nrow=nrow(mm)))
  }
  return(mm)
}


mztol <- function(mz,ppm,dmz){
  return(max(mz*ppm/1e6,dmz))
}


####Try to infer the position of the colun on the dataset
getHeader <- function(pt,mzname=c('mz','m.z.','MZ','mass','m/z'),
                      rtname=c('retention.time','rt','time'),
                      intname=c('int','intensity','area','height')){
  if(length(mzname)==0){
    mzname <- c('mz','m.z.','MZ','m/z','mass')
  }
  
  if(length(rtname)==0){
    rtname <- c('retention.time','rt','time')
  }
  
  if(length(intname)==0){
    intname <- c('int','intensity','area','height')
  }
  
  cnames <- tolower(trimws(colnames(pt)))
  
  pmz <- match(mzname,cnames)
  pok <- which(!is.na(pmz))
  if(length(pok)!=0){
    if(length(pok)>1){
      warning("Multiples columns found for m/z ",cnames[pmz[pok]],
              " selected: ",pmz[pok[1]])
    }
    pmz <- pmz[pok[1]]
  }
  
  prt <- match(rtname,cnames)
  pok <- which(!is.na(prt))
  if(length(pok)!=0){
    if(length(pok)>1){
      warning("Multiples columns found for rt ",cnames[prt[pok]],
              " selected: ",prt[pok[1]])
    }
    prt <- prt[pok[1]]
  }
  
  pint <- match(intname,cnames)
  pok <- which(!is.na(pint))
  if(length(pok)!=0){
    if(length(prt)>1){
      warning("Multiples columns found for intensity ",cnames[pint[pok]],
              " selected: ",pint[pok[1]])
    }
    pint <- pint[pok[1]]
  }
  return(c(pmz,prt,pint))
}


fortifyPeaktable <- function(peaktable,lar,supp_data=character(0),removeDuplicate=TRUE){
  ###A duplicate is found when two column have the same mass and retention time 
  ids <- paste(sprintf("%0.3f",peaktable[,1]),"-",sprintf("%0.2f",peaktable[,2]))
  
  selv <- !duplicated(ids)
  
  vm <- match(supp_data,colnames(peaktable))
  vm <- vm[!is.na(vm)]
  vpar <- params(lar)
  return(peaktable[selv,c(vpar[["col_mz"]],vpar[["col_rt"]],vpar[["col_int"]],vm),drop=FALSE])
}

fortifyPeaktableIndex <- function(peaktable,lar,supp_data=character(0),removeDuplicate=TRUE){
  ###A duplicate is found when two column have the same mass and retention time 
  ids <- paste(sprintf("%0.3f",peaktable[,1]),"-",sprintf("%0.2f",peaktable[,2]))
  
  selv <- !duplicated(ids)
  
  vm <- match(supp_data,colnames(peaktable))
  vm <- vm[!is.na(vm)]
  vpar <- params(lar)
  return(list(peaktable[selv,c(vpar[["col_mz"]],vpar[["col_rt"]],vpar[["col_int"]],vm),drop=FALSE],
         which(selv)))
}

get_os <- function() {
  if (.Platform$OS.type == "windows") {
    return("win")
  } else if (Sys.info()["sysname"] == "Darwin") {
    return("mac")
  } else if (.Platform$OS.type == "unix") {
    return("unix")
  } else {
    stop("Unknown OS")
  }
}



generateID <- function(lid,ndigits=3){
  # ff <- paste("%0.",ndigits,"f",sep="")
  return(paste(sprintf("%.f",mz*10^ndigits),sprintf("%.f",rt*10^ndigits),sample.int(10000,size=1),sep="_"))
}

write_dgCMatrix_csv <- function(mat,
                                filename,
                                chunk_size = 1000,
                                replace_old="NONE",replace_new=NA,
                                append = FALSE) {
  
  #library(Matrix)
  #library(data.table)
  
  # Transpose so retrieval of "rows" is much faster
  mat <- Matrix::t(mat)
  
  # Row names
  row_names <- colnames(mat)
  
  # gene names are now columns
  col_names <- rownames(mat)
  
  n_row <- ncol(mat)
  n_col <- length(col_names)
  
  n_chunks <- floor(n_row/chunk_size)
  
  # Initial chunk
  chunk <- 1
  chunk_start <- 1 + chunk_size * (chunk - 1)
  chunk_end <- chunk_size * chunk
  if(chunk_end>ncol(mat)){
    chunk_end <- ncol(mat)
  }
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  if(replace_old!="NONE"){
    chunk_mat[chunk_mat==replace_old] <- replace_new
  }
  
  chunk_df <- as.data.frame(chunk_mat)
  
  data.table::fwrite(chunk_df, file = filename, append = append)
  if(chunk_end < (chunk_size * chunk)){
    return(NULL)
  }
  
  # chunkation over chunks
  if(n_chunks>=2){
    for(chunk in 2:n_chunks) {
      chunk_start <- 1 + chunk_size * (chunk - 1)
      chunk_end <- chunk_size * chunk
      chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
      if(replace_old!="NONE"){
        chunk_mat[chunk_mat==replace_old] <- replace_new
      }
      chunk_df <- as.data.frame(chunk_mat)
      data.table::fwrite(chunk_df, file = filename, append = T)
    }
  }
  
  # Remaining samples
  chunk_start <- (n_chunks*chunk_size + 1)
  chunk_end <- n_row
  if(chunk_start>chunk_end){
    return(NULL)
  }
  chunk_mat <- t(as.matrix(mat[,chunk_start:chunk_end]))
  if(replace_old!="NONE"){
    chunk_mat[chunk_mat==replace_old] <- replace_new
  }
  chunk_df <- as.data.frame(chunk_mat)
  data.table::fwrite(chunk_df, file = filename, append = T)
}

