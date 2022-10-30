###online version using the onlineLCMSaligner method
suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(DBI,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(onlineLCMSaligner,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(data.table,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(stringr,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))

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
sink(file=stdout())

args <- commandArgs(trailingOnly = TRUE)
PATH_DB<- args[1]
PATH_BLOCKS <- args[2]
PATH_ALIGNMENT <- args[3]
PATH_OUT_DATAMATRIX <- args[4]
VAL_INTENSITY <- args[5]
RTTOL <- as.numeric(args[6])
MZTOL <- as.numeric(args[7])
MZPPM <- as.numeric(args[8])
NUM_REF <- as.numeric(args[9])
ALPHA_RT <- as.numeric(args[10])
NUM_CORES <- as.numeric(args[11])
# OUTFIGURE <- args[12]

###Creating the parallel porcessing objects.
bpp <- NULL
if(get_os()=="win"){
  bpp <- SnowParam(workers = NUM_CORES)
}else{
  bpp <- MulticoreParam(workers = NUM_CORES)
}
if(NUM_CORES==1){
    bpp <- SerialParam()
}


if(!file.exists(PATH_OUT_DATAMATRIX)){
    dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
    all_peaktables <- dbGetQuery(dbb, "SELECT output_ms FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")[, 1]
    dbDisconnect(dbb)
    
    ###We check the number of file by batch. No more than 700000 peaks
    to_check <- sample(all_peaktables,min(length(all_peaktables),10))
    MAX_PEAKS <- 700000
    peaks_by_sample <- mean(sapply(to_check,function(x){
      tt <- read.table(x,sep=",",header=TRUE)
      return(nrow(tt))
    }))
    
    max_by_batch <- floor(MAX_PEAKS/peaks_by_sample)
    ###We check the columns of the peaktable.
    pt <- read.table(all_peaktables[1],header=TRUE,sep=",")
    cnames <- colnames(pt)

    sel_args <- !cnames %in% c("mz","rt","rt_min","rt_max","mz_min","mz_max")
    supp_args <- cnames[sel_args]


    ##mz and rt are always stored
      lam <- LCMSAlignerModelFromDirectoryByBatch(all_peaktables,
               path_model=PATH_ALIGNMENT,
               output=PATH_BLOCKS,save_interval=max_by_batch,
               num_file=20,num_peaks=NUM_REF,col_int=VAL_INTENSITY,reset = FALSE,
               ppm = MZPPM, dmz=MZTOL,rt = RTTOL,rt_dens=RTTOL/2,n_clusters=10,
               supp_data=supp_args,ransac_l1=ALPHA_RT,bpp=bpp,
               max_cor=RTTOL*3,by_batch=max_by_batch,clustering=FALSE)
      if(is.null(lam)){
        stop("Empty datamatrix generated.")
      }

    ###We always filter out the peaks detected only once.
    vexp <- suppressMessages(suppressWarnings(exportDataMatrix(lam,path=PATH_OUT_DATAMATRIX,quant_var = VAL_INTENSITY,subvariable=which(lam@peaks$num>=2),summary_vars=c("mz","rt","rt_cor",supp_args))))

    if(VAL_INTENSITY %in% c("area","intensity","height")){
      ###In all case we reorder the datamatrix to avoid the first 
      dm <- fread(PATH_OUT_DATAMATRIX,sep="\t")
      ocnames <- colnames(dm)
      quant_prefix <- paste(str_split(ocnames[length(ocnames)],fixed("_"))[[1]][1],"_",sep="")
      quant_cols <- which(startsWith(ocnames,quant_prefix))
      quant_cols <- quant_cols[quant_cols>20] # removes the first ones
      raw_names_db <- basename(all_peaktables)
      raw_names_dm <- str_sub(colnames(dm)[quant_cols],nchar(quant_prefix)+1)
      # match(raw_names_dm,raw_names_db)
      new_idx <- quant_cols[match(raw_names_dm,raw_names_db)]
      ocnames[new_idx] <- lapply(ocnames[quant_cols],str_replace,quant_prefix,'quant_')
      dm[,new_idx] <- dm[,..quant_cols]
      colnames(dm) <- unlist(ocnames)
      fwrite(dm,file = PATH_OUT_DATAMATRIX,sep="\t")
    }
    
    
}else{
  cat("Data matrix already exists. Alignment won't be performed: ",PATH_OUT_DATAMATRIX,"\n",file=stdout())
}
