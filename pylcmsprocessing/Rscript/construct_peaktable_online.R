###online version using the onlineLCMSaligner method

suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(DBI,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(onlineLCMSaligner,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(jsonlite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))


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




groupPeaksDensity <- function (peaks, sampleGroups, bw = 10, minFraction = 0.5, minSamples = 1,
                               binSize = 0.01, maxFeatures = 50, sleep = 0)
{
  if (missing(sampleGroups))
    stop("Parameter 'sampleGroups' is missing! This should be a vector of ",
         "length equal to the number of samples specifying the group ",
         "assignment of the samples.")
  if (missing(peaks))
    stop("Parameter 'peaks' is missing!")
  if (!is.matrix(peaks) & !is.data.frame(peaks)){
    stop("'peaks' has to be a 'matrix' or a 'data.frame'!")
  }

  ###Store the current peak assignment
  peak_assign <- rep(NA_real_,nrow(peaks))
  colt <- c("time","rt")
  vrt <- match(colt,colnames(peaks))
  vrt <- colnames(peaks)[vrt[!sapply(vrt,is.na)][1]]
  .reqCols <- c("mz", vrt, "sample")
  if (!all(.reqCols %in% colnames(peaks)))
    stop("Required columns ", paste0("'", .reqCols[!.reqCols %in%
                                                     colnames(peaks)], "'", collapse = ", "), " not found in 'peaks' parameter")
  sampleGroups <- as.character(sampleGroups)
  sampleGroupNames <- unique(sampleGroups)
  sampleGroupTable <- table(sampleGroups)
  nSampleGroups <- length(sampleGroupTable)
  if (max(peaks[, "sample"]) > length(sampleGroups))
    stop("Sample indices in 'peaks' are larger than there are sample",
         " groups specified with 'sampleGroups'!")
  peakOrder <- order(peaks[, "mz"])
  peaks <- peaks[peakOrder, .reqCols, drop = FALSE]
  rownames(peaks) <- NULL
  rtRange <- range(peaks[, vrt])
  mass <- seq(peaks[1, "mz"], peaks[nrow(peaks), "mz"] + binSize,
              by = binSize/2)
  masspos <- xcms:::findEqualGreaterM(peaks[, "mz"], mass)
  groupmat <- matrix(nrow = 512, ncol = 7 + nSampleGroups)
  groupindex <- vector("list", 512)
  densFrom <- rtRange[1] - 3 * bw
  densTo <- rtRange[2] + 3 * bw
  densN <- max(512, 2 * 2^(ceiling(log2(diff(rtRange)/(bw/2)))))
  endIdx <- 0
  num <- 0
  gcount <- integer(nSampleGroups)
  for (i in seq_len(length(mass) - 2)) {
    startIdx <- masspos[i]
    endIdx <- masspos[i + 2] - 1
    if (endIdx - startIdx < 0)
      next
    curMat <- peaks[startIdx:endIdx, , drop = FALSE]
    den <- density(curMat[, vrt], bw = bw, from = densFrom,
                   to = densTo, n = densN)
    maxden <- max(den$y)
    deny <- den$y
    snum <- 0
    while (deny[maxy <- which.max(deny)] > maxden/20 && snum <
           maxFeatures) {
      grange <- xcms:::descendMin(deny, maxy)
      deny[grange[1]:grange[2]] <- 0
      gidx <- which(curMat[, vrt] >= den$x[grange[1]] &
                      curMat[, vrt] <= den$x[grange[2]])
      tt <- table(sampleGroups[unique(curMat[gidx, "sample"])])
      if (!any(tt/sampleGroupTable[names(tt)] >= minFraction &
               tt >= minSamples))
        next
      snum <- snum + 1
      num <- num + 1
      if (num > nrow(groupmat)) {
        groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat),
                                           ncol = ncol(groupmat)))
        groupindex <- c(groupindex, vector("list", length(groupindex)))
      }
      gcount <- rep(0, length(sampleGroupNames))
      names(gcount) <- sampleGroupNames
      gcount[names(tt)] <- as.numeric(tt)
      groupmat[num, 1] <- median(curMat[gidx, "mz"])
      groupmat[num, 2:3] <- range(curMat[gidx, "mz"])
      groupmat[num, 4] <- median(curMat[gidx, vrt])
      groupmat[num, 5:6] <- range(curMat[gidx, vrt])
      groupmat[num, 7] <- length(gidx)
      groupmat[num, 7 + seq(along = gcount)] <- gcount
      groupindex[[num]] <- sort(peakOrder[(startIdx:endIdx)[gidx]])
      peak_assign[peakOrder[(startIdx:endIdx)[gidx]]] <- num
    }
  }
  message("OK")
  colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed",
                          "rtmin", "rtmax", "npeaks", sampleGroupNames)
  groupmat <- groupmat[seq_len(num), , drop = FALSE]
  groupindex <- groupindex[seq_len(num)]
  numsamp <- rowSums(groupmat[, (match("npeaks", colnames(groupmat)) +
                                   1):ncol(groupmat), drop = FALSE])
  uorder <- order(-numsamp, groupmat[, "npeaks"])
  uindex <- xcms:::rectUnique(groupmat[, c("mzmin", "mzmax", "rtmin",
                                           "rtmax"), drop = FALSE], uorder)
  return(list(featureDefinitions = groupmat[uindex, , drop = FALSE],
              peakIndex = groupindex[uindex], peakAssign = peak_assign))
}





makeDataMatrix <- function(samples,sampnames,path_output_dm,path_output_index,funGroup=groupPeaksDensity,
                           int = "intensity",gp="group",  bw = 5, mztol = 0.005){
  suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
  suppressWarnings(suppressMessages(library(xcms,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
  suppressWarnings(suppressMessages(library(stringr,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))

  ###Construct peak table from porcessing table
  constructPeakTable <-  function(samples){
    ###We read all the table

    ###we tranform processing into a list to be passed to mapply
    vr <- mapply(samples,1:length(samples),FUN=function(x,y,intensity){
      if(!file.exists(x)){
        message("File ",x,"not found.")
        return(data.frame(mz=numeric(0),rt=numeric(0),intensity=numeric(0),
                          peakwidth=numeric(0),SN=numeric(0),right_on_left_assymetry=numeric(0)))
      }

      #"mz","rt","height","intensity","rt_min","rt_max","mz_min","mz_max","SN","peakwidth","right_on_left_assymetry"
      ####We read all the peak table
      tab <-  read.table(x,sep=",",header=TRUE)
      pint <- match(intensity,colnames(tab))

      if(is.na(pint)) message("intensity not found in file : ",x)
      ### The first two columns are always the mz and the rt
      tab <- tab[,c("mz","rt",intensity,"peakwidth","SN","right_on_left_assymetry"),drop=FALSE]
      tab$sample <-  rep(y,nrow(tab))
      return(tab)
    },SIMPLIFY = FALSE,MoreArgs=list("intensity"=int))

    ####Making the sample table.
    vr <-  do.call(rbind,vr)

    ####We remove all the peaks with NA in mass (may open for some open MS version)
    vr <- vr[!is.na(vr[,"mz"]),]

    return(vr)
  }

  message("processing ",basename(path_output_dm))
  ###We get the replicate of each sample
  pt <- constructPeakTable(samples)

  suppressMessages(mat <- funGroup(pt, rep(1,length(samples)), bw = bw, minFraction = 0.0, minSamples = 1,
                                   binSize = mztol, maxFeatures = 20, sleep = 0))
  message("Alignment finished!")
  ###We complete the data matrix using the pt table.
  maxSample <- length(samples)
  ###We build all the intensity table
  intTable <-  sapply(mat$peakIndex,function(x,pt,maxSample,missingValue,intensity){
    vline <- rep(missingValue,maxSample)
    vline[pt[x,"sample"]] <- pt[x,intensity]
    vline
  },pt=pt,maxSample=maxSample,missingValue=NA_real_,intensity=int,simplify=FALSE)

  ###We rewrite all the peaktables while storing the informations
  intTable <-  do.call(rbind,intTable)

  mean_peakwidth <- sapply(mat$peakIndex,function(x,pt){
    mean(pt[x,"peakwidth"],na.rm=TRUE)
  },pt=pt,simplify=FALSE)

  sd_peakwidth <- sapply(mat$peakIndex,function(x,pt){
    sd(pt[x,"peakwidth"],na.rm=TRUE)
  },pt=pt,simplify=FALSE)


  mean_assymetry <- sapply(mat$peakIndex,function(x,pt){
    median(pt[x,"right_on_left_assymetry"])
  },pt=pt,simplify=FALSE)

  sd_assymetry <- sapply(mat$peakIndex,function(x,pt){
    sd(pt[x,"right_on_left_assymetry"])
  },pt=pt,simplify=FALSE)


  mean_SN <- sapply(mat$peakIndex,function(x,pt){
    mean(pt[x,"SN"],na.rm=TRUE)
  },pt=pt,simplify=FALSE)

  ####Building the datamatrix using all the available information
  cnames <- c("mz","rt","mz_min","mz_max","rt_min","rt_max",'mean_peakwidth',
              'sd_peakwidth','median_assymetry','sd_assymetry',"mean_SN",paste("intensity",sampnames,sep="_"))
  dmm <- cbind(mat$featureDefinitions[,c("mzmed","rtmed","mzmin","mzmax","rtmin","rtmax")],
               mean_peakwidth,sd_peakwidth,mean_assymetry,sd_assymetry,mean_SN,intTable)
  colnames(dmm) <- cnames

  ###Writing the datamatrx
  write.table(dmm,file = path_output_dm,sep=",",row.names = FALSE)

  ###Writing the indices
  ff <- file(path_output_index,"w")
  write_json(mat$peakIndex,path = ff)
  close.connection(ff)
}

args <- commandArgs(trailingOnly = TRUE)

PATH_PEAKTABLES <- args[1]
PATH_BLOCKS <- args[2]
PATH_ALIGNMENT <- args[3]
PATH_OUT_DATAMATRIX <- args[4] #PATH_DATAMATRIX <- "test_datamatrix.csv"
PATH_OUT_IDX <- args[5] #PATH_DATAMATRIX <- "test_output.csv"
VAL_INTENSITY <- args[6]
RTTOL <- as.numeric(args[7])
MZTOL <- as.numeric(args[8])
MZPPM <- as.numeric(args[9])
NUM_REF <- as.numeric(args[10])
NUM_CLUSTERS <- as.numeric(args[11])


message("Beginning grouping using metric, ",VAL_INTENSITY)

lam <- LCMSAlignerModelFromDirectory(PATH_PEAKTABLES,
                      path_model=PATH_ALIGNMENT,
                       output=PATH_BLOCKS,save_interval=15,
                       num_file=20,num_peaks=NUM_REF,col_int=VAL_INTENSITY,reset = FALSE,allowed_jump=10,
                       ppm = MZPPM, dmz=MZTOL,rt = RTTOL,n_clusters=NUM_CLUSTERS)
##


bioldDataMatrix()

vfiles <- read.table(FNAMES,sep=",",header=FALSE,stringsAsFactors = FALSE)

tres <- makeDataMatrix(vfiles[,1],vfiles[,2],
                       PATH_DATAMATRIX,PATH_IDX,funGroup=groupPeaksDensity,int = VAL_INTENSITY,bw = bw, mztol=binsize)
message("Alignment done")
