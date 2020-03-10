library(BiocParallel,quietly=TRUE,warn.conflicts = FALSE)
library(MZmineXMLManipulator,quietly=TRUE,warn.conflicts = FALSE)
library(tools,quietly=TRUE,warn.conflicts = FALSE)

args <- commandArgs(trailingOnly=TRUE)
infiles <-  args[1]
MAX_WORKERS <- as.numeric(args[2])

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

bpp <- NULL
if(get_os()=="win"){
  bpp <- SnowParam(workers = MAX_WORKERS)
}else{
  bpp <- MulticoreParam(workers = MAX_WORKERS)
}

####Name of the file to convert
fnames <- readLines(infiles)

###The output name of the file is the same as the first one.
convPeakTableFromPath <-  function(pt){
  library(MZmineXMLManipulator,quietly=TRUE,warn.conflicts = FALSE)
  library(stringr,quietly=TRUE,warn.conflicts = FALSE)

  peaktable <- NULL
  peaktable <-  read.table(pt,header=TRUE,sep=",",stringsAsFactors = FALSE)
  peaktable <- tryCatch(convertPeakTable(peaktable),error=function(e){return(NA)})

  if(is.na(peaktable)) return(NULL)
  ###The file is written in place.
  write.table(peaktable,file = pt,row.names = FALSE,sep = ",")
  return(NULL)
}
bres <- bplapply(X=fnames,FUN=convPeakTableFromPath,BPPARAM = bpp)
