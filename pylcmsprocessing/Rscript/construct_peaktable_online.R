###online version using the onlineLCMSaligner method

options(error=traceback)

suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(DBI,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(onlineLCMSaligner,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))


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

args <- commandArgs(trailingOnly = TRUE)

PATH_DB<- args[1]
PATH_BLOCKS <- args[2]
PATH_ALIGNMENT <- args[3]
PATH_OUT_DATAMATRIX <- args[4] #PATH_DATAMATRIX <- "test_datamatrix.csv"
VAL_INTENSITY <- args[5]
RTTOL <- as.numeric(args[6])
MZTOL <- as.numeric(args[7])
MZPPM <- as.numeric(args[8])
NUM_REF <- as.numeric(args[9])
ALPHA_RT <- as.numeric(args[10])
NUM_CORES <- as.numeric(args[11])
OUTFIGURE <- args[12]

###Creating the parallel porcessing objects.
bpp <- NULL
if(get_os()=="win"){
  bpp <- SnowParam(workers = NUM_CORES)
}else{
  bpp <- MulticoreParam(workers = NUM_CORES)
}
if(!file.exists(PATH_OUT_DATAMATRIX)){
message("Beginning grouping using metric, ",VAL_INTENSITY)

dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
    all_peaktables <- dbGetQuery(dbb, "SELECT output_ms FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")[, 1]
dbDisconnect(dbb)

##mz and rt are always stored
lam <- LCMSAlignerModelFromDirectory(all_peaktables,
                      path_model=PATH_ALIGNMENT,
                       output=PATH_BLOCKS,save_interval=15,
                       num_file=20,num_peaks=NUM_REF,col_int=VAL_INTENSITY,reset = FALSE,
                       ppm = MZPPM, dmz=MZTOL,rt = RTTOL,rt_dens=RTTOL/2,n_clusters=10,
                       supp_data=c("peakwidth","SN","right_on_left_assymetry","height","intensity"),ransac_l1=ALPHA_RT,
                      max_cor=RTTOL*3,clustering=TRUE)

###We always remove single peaks.
if(!file.exists(OUTFIGURE)){
  pdf(OUTFIGURE)
  plotDevRt(lam,int_threshold=0.15)
  dev.off()
}else{
  message("Rt deviation figure already exists.")
}

###We always filter out the peaks detected only once.

  vexp <- exportDataMatrix(lam,path=PATH_OUT_DATAMATRIX,quant_var = VAL_INTENSITY,subvariable=which(lam@peaks$num>=2),summary_vars=c("mz","rt","rt_cor","peakwidth","SN","right_on_left_assymetry"))
  message("Alignment done")
}else{
  message("Data matrix already exists alignement won t be performed: ",PATH_OUT_DATAMATRIX)
}
