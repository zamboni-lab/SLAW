library(DBI)
library(RSQLite)
library(stringr)
library(openssl)
library(BiocParallel)

# INNER JOIN SYNTAX
# SELECT path,output_ms FROM samples INNER JOIN processing on sample = samples.id;

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



readMS2mzML <- function(path, outfile, noise_level=0) {
  suppressMessages(suppressWarnings(library(MSnbase, quietly = TRUE)))
  
  convertToMgfTxt <- function(raw_info, raw_spec, vid, noise_level=0) {
    if(nrow(raw_spec)==0) return(character(0))
    sel_peaks <- which(as.numeric(raw_spec[, 2])>noise_level)
    if(length(sel_peaks)==0) return(character(0))
    raw_spec <- raw_spec[sel_peaks,,drop=FALSE]
    raw_spec[, 1] <- sprintf("%0.4f", raw_spec[, 1])
    raw_spec[, 2] <- sprintf("%.3g", as.numeric(raw_spec[, 2]))
    txt_spec <- apply(raw_spec, 1, paste, collapse = " ")
    c(
      "BEGIN IONS",
      paste("FEATURE ID=", vid, sep = ""),
      paste("PEPMASS=", sprintf("%0.4f", raw_info$precursorMZ), sep = ""),
      paste("SCANS=", raw_info$precursorScanNum, sep = ""),
      paste(
        "RTINSECONDS=",
        sprintf("%0.2f", raw_info$retentionTime),
        sep = ""
      ),
      paste(
        "CHARGE=",
        raw_info$polarity,
        ifelse(sign(raw_info$polarity) == 1, "+", "-"),
        sep = ""
      ),
      paste(
        "COLLISIONENERGY=",
        raw_info$collisionEnergy,
        sep = ""
      ),
      "MSLEVEL=2",
      txt_spec,
      "END IONS"
    )
  }
  if(file.exists(outfile)) return(NA)
  
  of <- openMSfile(path)
  vhh <- header(of)
  rrh <- vhh[vhh$msLevel == 2, ]
  if(nrow(rrh)==0) return(NULL)
  all_peaks <- peaks(of, rrh$seqNum)
  if(!is.list(all_peaks)){
    all_peaks <- list(all_peaks)
  }
  tsplit <- split(rrh, f = 1:nrow(rrh))
  all_txt <-
    do.call(c, mapply(tsplit, all_peaks, as.list(seq_along(tsplit)), FUN =
                        convertToMgfTxt,MoreArgs=list(noise_level=noise_level),SIMPLIFY = FALSE))
  if(length(all_txt)==0) return(NULL)
  ff <- file(outfile)
  writeLines(all_txt, con = ff)
  close(ff)
}


args <- commandArgs(trailingOnly = TRUE)
# args <- c("U:/users/Alexis/data/slaw_evaluation/test_data/output/processing_db.sqlite",
#           "U:/users/Alexis/sandbox/test_output","500","4","TRUE")
PATH_DB <- args[1]
DIR_OUTPUT <- args[2]
NOISE_LEVEL <- as.numeric(args[3])
NUM_CORES <- as.integer(args[4])
ALL <- as.logical(args[5])

bpp <- NULL

if (get_os() == "win") {
  bpp <- SnowParam(workers = NUM_CORES, progressbar = TRUE)
} else{
  bpp <- MulticoreParam(workers = min(NUM_CORES, 4), progressbar = TRUE)
}



##SELECT ALL THE MS2FILE
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
infos <-
  dbGetQuery(dbb, "SELECT path,id FROM samples WHERE level='MS2'")
count_ms1 <-
  unlist(dbGetQuery(dbb, "SELECT COUNT(*) FROM samples WHERE level='MS1'"))
dbDisconnect(dbb)

raw_files <- infos[, 1]
# raw_files <- str_replace(raw_files,"/sauer1/","U:/")
raw_files <- c("U:/users/Alexis/examples_lcms_workflow/input_ms2/BEH30_2min_LipidMix_DDA_1.mzML",
               "U:/users/Alexis/examples_lcms_workflow/input_ms2/BEH30_2min_LipidMix_DDA_2.mzML")
vids <- as.integer(infos[, 2])
num_files <- length(raw_files)

need_computing <- numeric(0)

###We create all theoutput path
if(length(raw_files)>0){
all_path <-
  paste("MS2_",
        str_split(basename(raw_files), fixed("."), simplify = TRUE)[, 1],
        ".mgf",
        sep = "")
all_path <- file.path(DIR_OUTPUT, all_path)
all_hash <- md5(raw_files)

need_computing <- which(!file.exists(all_path))
}

if (length(need_computing) != 0) {
  vv <- bpmapply(
    raw_files[need_computing],
    all_path[need_computing],
    FUN = readMS2mzML,
    SIMPLIFY = FALSE,
    BPPARAM = bpp,
    MoreArgs=list(noise_level=NOISE_LEVEL)
  )
  
  ###We now add the output path to the processing table.
  supp_processing <-
    data.frame(
      id = as.integer(seq_along(raw_files) + count_ms1),
      peakpicking = as.integer(rep(1, num_files)),
      sample = vids,
      input = raw_files,
      hash = as.character(all_hash),
      output_ms = rep("NOT PROCESSED", num_files),
      output_ms2 = all_path,
      valid = as.integer(rep(1, num_files)),
      step = as.integer(rep(1, num_files))
    )
  
  supp_processing <- supp_processing[need_computing, , drop = FALSE]
  
  
  ####We now update the database.
  dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
  vev <- dbWriteTable(conn = dbb,
                       supp_processing,
                       name = "processing",append=TRUE)
  dbDisconnect(dbb)
}


if(ALL){
  message("Extracting all MS-MS spectra.")
  dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
  all_msms <- dbGetQuery(dbb, "SELECT path,output_ms2 FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1'")
  dbDisconnect(dbb)
  pout <- all_msms[,2]
  psamples <- all_msms[,1]
  need_computing <- which(!file.exists(pout))
  if (length(need_computing) != 0) {
    vv <- bpmapply(
      psamples[need_computing],
      pout[need_computing],
      FUN = readMS2mzML,
      SIMPLIFY = FALSE,
      BPPARAM = bpp,
      MoreArgs=list(noise_level=NOISE_LEVEL)
    )
  }
  message("MS-MS spectra extraction finished")
}

# readMS2mzML("U:/users/Alexis/data/slaw_evaluation/MTBLS865/mzML/neg_AEG15_Ac.mi._SPIL2_1-B,6_01_18851.mzML","U:/users/Alexis/data/slaw_evaluation/MTBLS865/neg_AEG15_Ac.mi._SPIL2_1-B,6_01_18851.mgf",noise_level=1000)
