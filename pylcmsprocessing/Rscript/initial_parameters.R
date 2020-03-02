suppressWarnings(suppressMessages(library(xcms,quietly=TRUE)))
suppressWarnings(suppressMessages(library(RSQLite,quietly=TRUE)))
suppressWarnings(suppressMessages(library(DBI,quietly=TRUE)))
suppressWarnings(suppressMessages(library(BiocParallel,quietly=TRUE)))
suppressWarnings(suppressMessages(library(yaml,quietly=TRUE)))
###Given a raw files, given a raw file, we estiamte the main parameters,
### The polarity and the rest of the data


args <- commandArgs(trailingOnly = TRUE)

# args <-
#   c(
#     "U:/users/Alexis/examples_lcms_workflow/output/processing_db.sqlite",
#     "U:/users/Alexis/examples_lcms_workflow/output/parameters.txt",
#     "E:/out_parameters.txt",
#     "4"
#   )

##
PATH_DB <- args[1]
INITIAL_PAR <- args[2]
OUTPUT_PAR <- args[3]
NUM_CORES <- as.numeric(args[4])


##Fixed parameters
##The widening of the estimated peakwidth
WIDENING <- c(0.2, 0.1)
##The minimum authorized peakwidth in seconds
MIN_PEAKWIDTH <- 1
##The maximum authorized peakwidth in seconds
MAX_PEAKWIDTH <- 60


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


####We get a list of all the raw files.
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
raw_files <- dbGetQuery(dbb, "SELECT path FROM samples")[, 1]
dbDisconnect(dbb)

bpp <- NULL
if (get_os() == "win") {
  bpp <- SnowParam(workers = NUM_CORES)
} else{
  bpp <- MulticoreParam(workers = min(NUM_CORES, 4))
}


###We try to estimate the noise

findPeakswidthParameters <-  function(xraw, widening, min_pw, max_pw){
  
  suppressWarnings(suppressMessages(library(xcms,quietly=TRUE)))
  suppressWarnings(suppressMessages(library(pracma,quietly=TRUE)))
  if(is.character(xraw)){
    xraw <-   suppressWarnings(suppressMessages(xcmsRaw(xraw)))
  }
  mzrange <- xraw@mzrange
  
  
  ###We extract all the 1 size Windows
  cut_mzs <- seq(mzrange[1] - 0.5, mzrange[2] + 0.5, by = 1)
  
  eic_mzs <- vector(mode = "list", length = length(cut_mzs) - 1)
  for (i in seq_along(eic_mzs)) {
    eic_mzs[[i]] <- cut_mzs[c(i, i + 1)]
  }
  
  
  extractBiggestPeak <- function(xraw, mzlim) {

    ttime <- xraw@scantime
    rre <- rawEIC(xraw, mzrange = mzlim)
    int_seq <- rre[[2]]
    if (all(int_seq == 0))
      return(c(NA, c(NA, NA), list()))
    # library(pracma)
    
    val_filter <- ceiling(length(int_seq) * 3 / 100)
    val_filter <- val_filter + (val_filter + 1) %% 2
    
    smoothed_seq <- savgol(int_seq, val_filter)
    # plot(smoothed_seq,type="l")
    
    ##We detect the biggest peak
    
    pmm <- which.max(smoothed_seq)
    pmin <- pmm
    pmax <- pmm
    while (pmin > 1 &&
           smoothed_seq[pmin - 1] < smoothed_seq[pmin])
      pmin <- pmin - 1
    while (pmax < length(smoothed_seq) &&
           smoothed_seq[pmax + 1] < smoothed_seq[pmax])
      pmax <- pmax + 1
    
    if (pmax == length(int_seq))
      return(c(NA, c(NA, NA), list()))
    dseq <- diff(int_seq)
    ###We adjust the minimum and the maximum to find the first point where 2 points are increasing
    while (pmin < pmm && !(dseq[pmin] > 0.1 & dseq[pmin + 1] > 0.1)) {
      pmin <- pmin + 1
    }
    while (pmax > pmm && !(dseq[pmax + 1] < (-0.1) &
                           dseq[pmax] < (-0.1))) {
      pmax <- pmax - 1
    }
    pmin <- pmin
    pmax <- pmax
    return(list(int_seq[pmm], ttime[c(pmin, pmax)], c(pmin, pmax), int_seq[pmin:pmax]))
  }
  
  vee <-
    lapply(eic_mzs,
             FUN = extractBiggestPeak,
             xraw = xraw)
  vsel <- sapply(vee, function(x) {
    !is.na(x[1]) && (diff(x[[2]]) > 0.8)
  })
  vee <- vee[vsel]
  
  
  ###We get the peakwidth
  vpeakwidth <- sapply(vee, function(x) {
    x[[2]][2] - x[[2]][1]
  })
  
  
  ###The bandidth is always 0.2
  vden <- density(vpeakwidth, bw = 0.2)
  
  pmm <- which.max(vden$y)
  pmin <- pmax <- pmm
  while (pmin > 2 && vden$y[pmin - 1] < vden$y[pmin])
    pmin <- pmin - 1
  while (pmax < (length(vden$y) - 1) &&
         vden$y[pmax] > vden$y[pmax + 1])
    pmax <- pmax + 1

  
  ###We tae the density peak as the number of peaks.
  peakwidth <- c(vden$x[pmin], vden$x[pmax])
  enlargement_left <- peakwidth[1] * (widening[1])
  enlargement_right <- diff(peakwidth) * widening[2]
  peakwidth <-
    c(peakwidth[1] - enlargement_left, peakwidth[2] + enlargement_right)
  peakwidth <-
    c(max(peakwidth[1], min_pw),
      min(peakwidth[2], max_pw))
  
  ###We compute the time spent on a scan
  time_scan <-
    mean(sapply(vee, function(x) {
      diff(x[[2]] / diff(x[[3]]))
    }))
  
  ###We convert this knowledge in number of scans.
  num_scans <- floor(peakwidth[1] / time_scan)
  
  polarity <- unique(as.character(xraw@polarity))[1]
  if(polarity=="positive"){
    polarity <- 1
  }else{
    polarity <- -1
  }
  return(c(peakwidth,num_scans,polarity))
}

handlers <- suppressWarnings(suppressMessages(lapply(raw_files,xcmsRaw)))
resVal <- bplapply(handlers,FUN=findPeakswidthParameters,widening=WIDENING,min_pw=MIN_PEAKWIDTH,max_pw=MAX_PEAKWIDTH,BPPARAM = bpp)
estimated_par <- apply(do.call(rbind,resVal),2,median)
polarity = estimated_par[3]
if(polarity==1){
  Sys.setenv("POLARITY"="positive")
}else{
  Sys.setenv("POLARITY"="negative")
}

##We read the initial parameter estimation
params <- yaml.load_file(INITIAL_PAR)
params$peakpicking$peaks_deconvolution$peak_width_min$value <- estimated_par[1]/60
params$peakpicking$peaks_deconvolution$peak_width_max$value <- estimated_par[2]/60
params$peakpicking$traces_construction$min_scan$value <- as.integer(estimated_par[3])
write_yaml(params,OUTPUT_PAR)



