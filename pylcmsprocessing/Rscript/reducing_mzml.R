suppressWarnings(suppressMessages(library(xcms,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(mzR,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(stringr,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))


sink(file=stdout())
###Trying to extract region of interest rom a chormatogram.
args <- commandArgs(trailingOnly = TRUE)
###We just replace the input file if the path output is not given


INPUT_DIR <- args[1]

if(length(args)>1){
  PATH_OUTPUT <- args[2]
}

MAX_ISO <- 4
MAX_SEL <- 20

all_files <- list.files(INPUT_DIR,full.names=TRUE,pattern="mzX?ML")


extract_mzbins <- function(x,MAX_ISO,MAX_SEL){
  suppressWarnings(suppressMessages(library(xcms,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
  xraw <- suppressWarnings(suppressMessages(xcmsRaw(x)))
  rmz <- range(xraw@env$mz)
  mz_min <- floor(rmz[1])-0.5
  mz_max <- ceiling(rmz[2])+0.5
  prof <- suppressWarnings(suppressMessages(profMat(xraw,mzrange=c(mz_min,mz_max))))
  int_val <- suppressWarnings(suppressMessages(apply(prof,1,max)))
  selected <- rep(FALSE,length(int_val))
  ###We just select the row awith max intensity and the +3  to +4
  num_sel <- 1
  MAX_ISO <-5
  while(num_sel<MAX_SEL){
    sel_idx <- which(!selected)
    pmax <- sel_idx[which.max(int_val[sel_idx])]
    max_iso <- min(pmax+MAX_ISO,length(selected))
    pisos <- pmax:max_iso
    selected[pisos] <- TRUE  
    num_sel <- num_sel+1
  }
  return(list(mz_min,mz_max,selected))
}


sel_data <- bplapply(as.list(all_files),extract_mzbins,MAX_ISO=MAX_ISO,MAX_SEL=MAX_SEL)

min_mz_bin <- min(sapply(sel_data,"[[",i=1))
max_mz_bin <- max(sapply(sel_data,"[[",i=2))

sel <- rep(FALSE,max_mz_bin-min_mz_bin+1)

for(idx in seq_along(sel_data)){
  vmin <- sel_data[[idx]][[1]]
  sel[vmin-min_mz_bin+which(sel_data[[idx]][[3]])] <- TRUE
}

selected <- which(sel)

cat("Conserving ",length(selected)," m/z bins for optimization.",file=stdout())

mzlim_data <- cbind(min_mz_bin+selected-1-0.5,min_mz_bin+selected-0.5)

filter_mzml <- function(praw,poutput,mzlims){
  suppressWarnings(suppressMessages(library(xcms,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
  suppressWarnings(suppressMessages(library(mzR,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
  filter_mz_interval <- function(x, mzlims, msLevel., updatePeaksCount = TRUE){
      if (!missing(msLevel.)) {
        if (!(msLevel(x) %in% msLevel.)) 
          return(x)
      }
      sel <- rep(FALSE,length(x@mz))
      
      current_idx <- 1
      cmzmin <- mzlims[current_idx,1]
      cmzmax <- mzlims[current_idx,2]
      
      for(idx in seq_along(x@mz)){
        cmz <- x@mz[idx]
        if(cmz<cmzmin) next
        if(cmz>cmzmax){
          current_idx <- current_idx+1
          if(current_idx>nrow(mzlims)) break
          cmzmin <- mzlims[current_idx,1]
          cmzmax <- mzlims[current_idx,2]
          next
        }
        if((cmz<=cmzmax)&(cmz>=cmzmin)){
          sel[idx] <- TRUE
        }
      }
      
      if (!length(sel)) {
        msg <- paste0("No data points between ", mzmin, 
                      " and ", mzmax, " for spectrum with acquisition number ", 
                      acquisitionNum(x), ". Returning empty spectrum.")
        warning(paste(strwrap(msg), collapse = "\n"))
        x@mz <- x@intensity <- numeric()
        x@tic <- integer()
        x@peaksCount <- 0L
        return(x)
      }
      x@mz <- x@mz[sel]
      x@intensity <- x@intensity[sel]
      if (updatePeaksCount) 
        x@peaksCount <- as.integer(length(x@intensity))
      return(x)
  }
  
  xrr <-readMSData(praw,mode="inMemory",msLevel. = 1)
  xrr@assayData <- suppressWarnings(suppressMessages(as.environment(eapply(assayData(xrr),filter_mz_interval,mzlims=mzlim_data,msLevel.=1))))
  if(file.exists(poutput))
    file.remove(poutput)
  suppressMessages(suppressWarnings(writeMSData(xrr,file = poutput,outformat="mzml")))
}

nn <- suppressMessages(suppressWarnings(mapply(all_files,all_files,FUN=filter_mzml,MoreArgs=list(mzlims=mzlim_data))))
