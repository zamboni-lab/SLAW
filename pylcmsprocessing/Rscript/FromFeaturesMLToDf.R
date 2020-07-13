#PATH_TEST <- "C:/Users/dalexis/Documents/LCMS_processing/dev/LCMS_script/features_ex.featureXML"
args <-  commandArgs(trailingOnly = TRUE)

# args <- c("U:/users/Alexis/examples_lcms_workflow/tt4.featureML","U:/users/Alexis/examples_lcms_workflow/tt4.csv")

suppressMessages(library(XML))
suppressMessages(library(stringr))
suppressMessages(library(DBI))
suppressMessages(library(RSQLite))
# 
# args <- c("U:/users/Alexis/data/slaw_evaluation/MTBLS865/output/OPENMS/peaktables/peaktable_527_cbfc76e654c8cdb843043e670899fede.csv","U:/users/Alexis/data/slaw_evaluation/MSV83010/output_sub/OPENMS/peaktables/tout.csv")


feature_ml <- args[1]
outfile <-  args[2]


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
# mz	rt	height	intensity	rt_min	rt_max	mz_min	mz_max	SN	peakwidth	right_on_left_assymetry

convertXMLopenMS <-  function(fv) {
  suppressMessages(library(XML))
  suppressMessages(library(stringr))
  
  
  rtree <- xmlParse(fv,asText = FALSE)
  rootnode <- xmlRoot(rtree)
  rootsize <- xmlSize(rootnode)
  
  ##Detemrining the position of the root
  posroot <- which(names(rootnode)=="featureList")

  num_features <- xmlSize(rootnode[[posroot]])

  convertValue <- function(x) {
    istable <- FALSE
    
    val <- str_remove_all(x, "[\\[|\\]]")
    if (nchar(val) != nchar(x)) {
      istable <- TRUE
    }
    if (istable) {
      val <- str_split(val, pattern = fixed(","), simplify = TRUE)
    }
    
    res <- suppressWarnings(as.numeric(val))
    if (any(is.na(res)))
      res <- val
    return(res)
  }

  rts <- sapply(1:num_features,function(ix,rn){
    x <- rn[[ix]]
    convertValue(xmlGetAttr(x[[12]], "value"))
  },rn=rootnode[[posroot]])
  rts <- sort(unique(unlist(rts)))
  
  margin <- 0.002
  
  dfl <- sapply(1:num_features, function(ix, rn,margin=0.002) {
    x <- rn[[ix]]
    
    intensity <- x[[3]]
    allint <- convertValue(xmlGetAttr(x[[11]], "value"))
    allrt <- convertValue(xmlGetAttr(x[[12]], "value"))
    allmz <- convertValue(xmlGetAttr(x[[13]], "value"))
    num_signals <- length(allint)
    fwhm <- convertValue(xmlGetAttr(x[[9]], "value"))
    rt_min <- allrt-fwhm*0.75
    rt_max <- allrt+fwhm*0.75
    mz_min <- allmz-margin
    mz_max <- allmz+margin
    quality <- as.numeric(xmlValue(x[[6]]))
    fwhm <- rep(fwhm,length(allrt))
    quality <- rep(quality,length(allrt))
    
    tempdf <- data.frame(mz = allmz,rt = allrt,height=allint,
                         intensity = allint,rt_min=rt_min,rt_max=rt_max,
                         mz_min=mz_min,mz_max=mz_max,peakwidth=fwhm,
                         quality=quality)
    return(tempdf)
  }, rn = rootnode[[posroot]], margin=margin, simplify = FALSE)
  
  # mz	rt	height	intensity	rt_min	rt_max	mz_min	mz_max	SN	peakwidth	right_on_left_assymetry
  
  dfl <- do.call(rbind, dfl)
  
  dfl <-  dfl[order(dfl$mz, decreasing = FALSE), ]
  
  ###We map the retention time to the TRUE retention time.
  map_max <- sapply(dfl$rt_max,function(x,rts){which.min(abs(x-rts))},rts=rts)
  map_max <- mapply(map_max,as.list(dfl$rt_max),FUN=function(x,rt,rts){
    if(rt>rts[x]&x<length(rt)){
      x <- x+1
    }
    x
  },MoreArgs = list(rts=rts))
  map_min <- sapply(dfl$rt_min,function(x,rts){which.min(abs(x-rts))},rts=rts)
  map_min <- mapply(map_min,as.list(dfl$rt_min),FUN=function(x,rt,rts){
    if(rt<rts[x]&x>1){
      x <- x-1
    }
    x
  },MoreArgs = list(rts=rts))
  
  dfl$rt_min <- rts[map_min]/60
  dfl$rt_max <- rts[map_max]/60
  dfl$rt <- dfl$rt/60
  dfl$peakwidth <- dfl$peakwidth/60
  ####We dont remove the ifle as we need ti trestart the calculations eventually.
  return(dfl)
}
###If an error occurs we just output an emtpy data frame
converted_table <- tryCatch(convertXMLopenMS(feature_ml),error=function(e){
  return(data.frame(mz = numeric(0),rt = numeric(0),height=numeric(0),
                    intensity = numeric(0),rt_min=numeric(0),rt_max=numeric(0),
                    mz_min=numeric(0),mz_max=numeric(0),peakwidth=numeric(0),
                    quality=numeric(0)))
})
wt <- file.remove(feature_ml)
wt <-write.table(converted_table, file = outfile, row.names = FALSE,sep=",",col.names = TRUE)
