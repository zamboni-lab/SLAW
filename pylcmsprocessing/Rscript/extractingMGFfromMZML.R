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
    sel_peaks <- which(as.numeric(raw_spec[, 2])>noise_level)
    if(length(sel_peaks)==0) return(character())
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
      "MSLEVEL=2",
      txt_spec,
      "END IONS"
    )
  }
  
  of <- openMSfile(path)
  vhh <- header(of)
  rrh <- vhh[vhh$msLevel == 2, ]
  all_peaks <- peaks(of, rrh$seqNum)
  tsplit <- split(rrh, f = 1:nrow(rrh))
  all_txt <-
    do.call(c, mapply(tsplit, all_peaks, as.list(seq_along(tsplit)), FUN =
                        convertToMgfTxt,MoreArgs=list(noise_level=noise_level)))
  ff <- file(outfile)
  writeLines(all_txt, con = ff)
  close(ff)
}


args <- commandArgs(trailingOnly = TRUE)
# args <- c("U:/users/Alexis/examples_lcms_workflow/output/processing_db.sqlite",
#           "U:/users/Alexis/examples_lcms_workflow/tout_ms2","4")
PATH_DB <- args[1]
DIR_OUTPUT <- args[2]
NOISE_LEVEL <- as.numeric(args[3])
NUM_CORES <- as.integer(args[4])

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
vids <- as.integer(infos[, 2])
num_files <- length(raw_files)

###We create all theoutput path
all_path <-
  paste("MS2_",
        str_split(basename(raw_files), fixed("."), simplify = TRUE)[, 1],
        ".mgf",
        sep = "")
all_path <- file.path(DIR_OUTPUT, all_path)
all_hash <- md5(raw_files)

need_computing <- which(!file.exists(all_path))

if (length(need_computing) != 0) {
  bpmapply(
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
