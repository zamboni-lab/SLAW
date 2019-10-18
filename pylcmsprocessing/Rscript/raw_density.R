
groupPeaks <- function (peaks, sampleGroups, bw = 10, minFraction = 0.5, minSamples = 1, 
          binSize = 0.01, maxFeatures = 50, sleep = 0) 
{
  if (missing(sampleGroups)) 
    stop("Parameter 'sampleGroups' is missing! This should be a vector of ", 
         "length equal to the number of samples specifying the group ", 
         "assignment of the samples.")
  if (missing(peaks)) 
    stop("Parameter 'peaks' is missing!")
  if (!is.matrix(peaks) | is.data.frame(peaks)) 
    stop("'peaks' has to be a 'matrix' or a 'data.frame'!")
  .reqCols <- c("mz", "rt", "sample")
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
  rtRange <- range(peaks[, "rt"])
  mass <- seq(peaks[1, "mz"], peaks[nrow(peaks), "mz"] + binSize, 
              by = binSize/2)
  masspos <- xcms <- findEqualGreaterM(peaks[, "mz"], mass)
  groupmat <- matrix(nrow = 512, ncol = 7 + nSampleGroups)
  groupindex <- vector("list", 512)
  densFrom <- rtRange[1] - 3 * bw
  densTo <- rtRange[2] + 3 * bw
  densN <- max(512, 2 * 2^(ceiling(log2(diff(rtRange)/(bw/2)))))
  endIdx <- 0
  num <- 0
  gcount <- integer(nSampleGroups)
  message("Processing ", length(mass) - 1, " mz slices ... ", 
          appendLF = FALSE)
  for (i in seq_len(length(mass) - 2)) {
    startIdx <- masspos[i]
    endIdx <- masspos[i + 2] - 1
    if (endIdx - startIdx < 0) 
      next
    curMat <- peaks[startIdx:endIdx, , drop = FALSE]
    den <- density(curMat[, "rt"], bw = bw, from = densFrom, 
                   to = densTo, n = densN)
    maxden <- max(den$y)
    deny <- den$y
    snum <- 0
    while (deny[maxy <- which.max(deny)] > maxden/20 && snum < 
           maxFeatures) {
      grange <- xcms:::descendMin(deny, maxy)
      deny[grange[1]:grange[2]] <- 0
      gidx <- which(curMat[, "rt"] >= den$x[grange[1]] & 
                      curMat[, "rt"] <= den$x[grange[2]])
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
      groupmat[num, 4] <- median(curMat[gidx, "rt"])
      groupmat[num, 5:6] <- range(curMat[gidx, "rt"])
      groupmat[num, 7] <- length(gidx)
      groupmat[num, 7 + seq(along = gcount)] <- gcount
      groupindex[[num]] <- sort(peakOrder[(startIdx:endIdx)[gidx]])
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
              peakIndex = groupindex[uindex]))
}