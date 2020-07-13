suppressWarnings(suppressMessages(library(xcms,quietly=TRUE)))
suppressWarnings(suppressMessages(library(RSQLite,quietly=TRUE)))
suppressWarnings(suppressMessages(library(DBI,quietly=TRUE)))
suppressWarnings(suppressMessages(library(BiocParallel,quietly=TRUE)))
suppressWarnings(suppressMessages(library(yaml,quietly=TRUE)))
###Given a raw files, given a raw file, we estiamte the main parameters,
### The polarity and the rest of the data


args <- commandArgs(trailingOnly = TRUE)
# args <- c("U:/users/Alexis/data/slaw_evaluation/test_data/output/processing_db.sqlite",
#           "U:/users/Alexis/data/slaw_evaluation/test_data/output/processing_db.sqlite")
# 
# 
# ##
PATH_DB <- args[1]
INITIAL_PAR <- args[2]
OUTPUT_PAR <- args[3]
NUM_CORES <- as.numeric(args[4])
num_sel <- max(10,2*NUM_CORES)
message("Initial parameters estimation")

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
raw_files <- dbGetQuery(dbb, "SELECT path FROM samples WHERE types='QC'")[, 1]
if(length(raw_files)==0) raw_files <- dbGetQuery(dbb, "SELECT path FROM samples WHERE types='sample'")[, 1]
dbDisconnect(dbb)



# if(FALSE){
#   raw_files <- list.files("U:/users/Alexis/data/slaw_evaluation/optim_MTBLS1129/mzML",full.names = TRUE)
#   INITIAL_PAR <- "U:/users/Alexis/data/slaw_evaluation/optim_MTBLS1129/output/parameters.txt"
#   OUTPUT_PAR <- "U:/users/Alexis/data/slaw_evaluation/optim_MTBLS1129/output/out_parameters.txt"
#   NUM_CORES <- 4
#   num_sel <- length(raw_files)  
# }

if(num_sel>length(raw_files)){
  num_sel <- length(raw_files)
}
raw_files <- sample(raw_files,num_sel)

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

  ###We get the most intense scan
  best_time <- xraw@scantime[which.max(xraw@tic)]
  
  extractBiggestPeak <- function(xraw, mzlim) {

    ttime <- xraw@scantime
    rre <- rawEIC(xraw, mzrange = mzlim)
    int_seq <- rre[[2]]
    if (all(int_seq == 0))
      return(c(NA, c(NA, NA), list()))
    # library(pracma)
    
    val_filter <- ceiling(length(int_seq) * 1 / 100)
    val_filter <- val_filter + (val_filter + 1) %% 2
    smoothed_seq <- savgol(int_seq, val_filter)

    ##We detect the biggest peak
    pmm <- which.max(smoothed_seq)
    pmin <- pmm
    pmax <- pmm
    while (pmin > 1 &&
           int_seq[pmin - 1] < int_seq[pmin])
      pmin <- pmin - 1
    while (pmax < length(int_seq) &&
           int_seq[pmax + 1] < int_seq[pmax])
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
  num_scans <- max(ceiling(peakwidth[1] / time_scan),2)
  
  polarity <- unique(as.character(xraw@polarity))[1]
  if(polarity=="positive"){
    polarity <- 1
  }else{
    polarity <- -1
  }
  min_int <- min(xraw@env$intensity)
  max_time <- max(xraw@scantime)
  return(c(peakwidth,num_scans,polarity,min_int,best_time,max_time))
}

handlers <- suppressWarnings(suppressMessages(lapply(raw_files,xcmsRaw)))
# r1 <- findPeakswidthParameters(handlers[[1]],widening=WIDENING,min_pw=MIN_PEAKWIDTH,max_pw=MAX_PEAKWIDTH)

resVal <- bplapply(handlers,FUN=findPeakswidthParameters,widening=WIDENING,min_pw=MIN_PEAKWIDTH,max_pw=MAX_PEAKWIDTH,BPPARAM = bpp)
estimated_params <- do.call(rbind,resVal)
estimated_par <- apply(estimated_params,2,median)
polarity <- estimated_par[4]
total_comb <- sum(1:(length(handlers)-1))
frac <- sum(1:(ceiling(length(handlers)/2-1)))/total_comb

maxRtDev <- 2*quantile(abs(combn(estimated_params[,6],2,diff)),probs = frac)/60
max_time <- estimated_par[7]/100
maxRtDev <- max(max_time/60,maxRtDev)

##We read the initial parameter estimation
params <- yaml.load_file(INITIAL_PAR)
bmin <- estimated_par[1]/60
bmax <- estimated_par[2]/60

min_scan <- max(floor(estimated_par[3]),2)


###We check if it is TOF or Orbitrap
noiseMS1 <- max(estimated_par[5]*1.5,200)
noiseMS2 <- max(estimated_par[5]*1.5,200)
feature_limits <- noiseMS1*10
if(estimated_par[5]<500){
  ppm <- 7
  ppm_range <- c(3,40)
  dmz <- 0.003
  dmz_range <- c(0.002,0.008)
  bin_group <- c(0.004)
  bin_group_range <- c(0.002,0.01)
  feature_limits <- 2000
}else{  ###Orbitrap
  ppm <- 3
  ppm_range <- c(3,25)
  dmz <- 0.002
  dmz_range <- c(0.001,0.005)
  bin_group <- c(0.004)
  bin_group_range <- c(0.002,0.01)
  feature_limits <- 150000
}


params$peakpicking$peaks_deconvolution$peak_width$value <- c(bmin,bmax)
bbmin <- bmin*0.5
if(bbmin<0.01) bbmin <- 0.01

###peakpicking parameters
if("range" %in% names(params$peakpicking$peaks_deconvolution$peak_width)){
  params$peakpicking$peaks_deconvolution$peak_width$range$min <- c(bbmin,bmin*c(3))
  params$peakpicking$peaks_deconvolution$peak_width$range$max <- bmax*c(0.5,3)
}
if("range" %in% names(params$peakpicking$traces_construction$ppm)){
  params$peakpicking$traces_construction$ppm$value <- ppm
  params$peakpicking$traces_construction$ppm$range <- ppm_range
}
if("range" %in% names(params$peakpicking$traces_construction$dmz)){
  params$peakpicking$traces_construction$dmz$value <- dmz
  params$peakpicking$traces_construction$dmz$range <- dmz_range
}
if("range" %in% names(params$peakpicking$traces_construction$min_scan)){
  params$peakpicking$traces_construction$min_scan$value <- min_scan
  temp_range <- ceiling(min_scan*c(0.75,1))
  if(temp_range[1]==temp_range[2]){
    temp_range[2] <- temp_range[2]+2
  }
  params$peakpicking$traces_construction$min_scan$range <- temp_range
  
  if("range" %in% names(params$peakpicking$traces_construction$num_outliers)){
    params$peakpicking$traces_construction$num_outliers$value <- ceiling(0.75*min_scan)
    temp_range <- ceiling(min_scan*c(0.75,1))
    if(temp_range[1]==temp_range[2]){
      temp_range[2] <- temp_range[2]+2
    }
    params$peakpicking$traces_construction$num_outliers$range <- temp_range
  }
}


###Noise level
params$peakpicking$noise_level_ms1$value <- noiseMS1
params$peakpicking$noise_level_ms2$value <- noiseMS2

params$peakpicking$peaks_deconvolution$noise_level$value <- feature_limits
###Grouping parameters
if("range" %in% names(params$grouping$drt)){
  ###Max RT dev should be 0.5% of gradient length in every cases
  params$grouping$drt$range <- maxRtDev*c(0.1,0.4)
  params$grouping$drt$value <- maxRtDev*c(0.2)
}
  
  
if("range" %in% names(params$grouping$dmz)){
  ###Max RT dev should be 0.5% of gradient length in every cases
  params$grouping$dmz$range <- bin_group
  params$grouping$dmz$value <- bin_group_range
}

if("range" %in% names(params$grouping$ppm)){
  params$grouping$ppm$range <- ppm_range
  params$grouping$ppm$value <- ppm
}

if("range" %in% names(params$grouping$dmz)){
  params$grouping$dmz$range <- dmz_range
  params$grouping$dmz$value <- dmz
}


write_yaml(params,OUTPUT_PAR)



