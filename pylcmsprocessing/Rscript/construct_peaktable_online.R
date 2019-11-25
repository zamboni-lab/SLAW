###online version using the onlineLCMSaligner method

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
message("args: ",paste(args,collapse="_"))


PATH_PEAKTABLES <- args[1]
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

message("Beginning grouping using metric, ",VAL_INTENSITY)

##mz and rt are always stored
lam <- LCMSAlignerModelFromDirectory(PATH_PEAKTABLES,
                      path_model=PATH_ALIGNMENT,
                       output=PATH_BLOCKS,save_interval=15,
                       num_file=20,num_peaks=NUM_REF,col_int=VAL_INTENSITY,reset = FALSE,
                       ppm = MZPPM, dmz=MZTOL,rt = RTTOL,rt_dens=RTTOL/2,n_clusters=10,
                       supp_data=c("peakwidth","SN","right_on_left_assymetry"),ransac_l1=ALPHA_RT,
                      max_cor=RTTOL*3,clustering=FALSE)

###We always remove single peaks.
pdf(OUTFIGURE)
plotDevRt(lam,int_threshold=0.15)
dev.off()

###We always filter the peaks present a a single time.
if(!file.exists(PATH_OUT_DATAMATRIX)){
  r_datamatrix <- buildDataMatrix(lam,subvariable=which(lam@peaks$num>=2),summary_vars=c("mz","rt","rt_cor","peakwidth","SN","right_on_left_assymetry"))
  write.table(r_datamatrix,file = PATH_OUT_DATAMATRIX,sep=",",row.names = FALSE)
message("Alignment done")
}
