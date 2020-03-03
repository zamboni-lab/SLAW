#' @include allGenerics.R
#' @include classes.R

#' Creation of an LCMSAlignerFileHandler object
#'
#' @return An empty LCMSAlignerFileHandler object 
#' @export
#'
#' @examples
#' print('Examples to be put here')
LCMSAlignerFileHandler <- function(){
  new_env <- new.env(hash = TRUE, parent = parent.frame(), size = 100L)
  inv_env <- new.env(hash = TRUE, parent = parent.frame(), size = 100L)
  order_file <- new.env(hash = TRUE, parent = parent.frame(), size = 100L)
  temp <- new("LCMSAlignerFileHandler")
  temp@files <- new_env
  temp@inv_files <- inv_env
  temp@id_table <- order_file
  temp@num_files <- 0
  return(temp)
}

hashFile <- function(x){
  return(md5sum(x))
}


#' Add a file to LCMSALignerFIleHandler
#'
#' @param lafh An LCMSAligenrFileHandler object
#' @param path The path to the file ot be added
#'
#' @return  The filled LCMSFIleHAandler object
#' @export
#'
#' @examples
#' print('Examples ot be put here')
setMethod("addFile","LCMSAlignerFileHandler",function(object,path) {
  if(!file.exists(path)){
    stop('File does not exists, you can only add existing file')
  }
  vhash <- hashFile(path)
  val <- object@files[[vhash]]
  if(!is.null(val)){
    warning('File ',path,' has already been processed.')
    return(list(obj=object,added=FALSE, id=NA_real_))
  }else{
    object@num_files <- object@num_files+1
    assign(vhash,path,envir = object@files)
    assign(path,vhash,envir = object@inv_files)
    assign(vhash,object@num_files,envir = object@id_table)
    return(list(obj=object,added=TRUE, id=object@num_files))    
  }
})

assignFile <- function(venv,path,vmd5=NULL){
  if(!file.exists(path)){
    stop('File does not exists')
  }
  if(is.null(vmd5)) vmd5 <- md5sum(path)
  assign(vmd5,value = path,venv)
}

setMethod("length","LCMSAlignerFileHandler",function(x){
  return(length(x@files))
})

###Return the numeric id of a file form his path
getId <- function(lafh,path){
  vhash <- hashFile(path)
  return(lafh@id_table[[vhash]]) 
}


setMethod("isProcessed","LCMSAlignerFileHandler",function(object,paths){
  info_pro <- logical(length(paths))
  for(ip in seq_along(paths)){
    info_pro[ip] <- !is.null(object@inv_files[[paths[ip]]])
  }
  return(info_pro)
})

##Return the basename of the raw file in the correct order of id
getAllFiles <- function(lafh){
  vids <- unname(unlist(ls(lafh@files)))
  vfiles <- sapply(vids,function(x,envv){envv[[x]]},envv=lafh@files)
  vpos <- sapply(vids,function(x,envv){envv[[x]]},envv=lafh@id_table)
  cvfiles <- character(length(vfiles))
  cvfiles[vpos] <- vfiles
  return(basename(cvfiles))
  
}