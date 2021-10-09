###Gap filling using an alignment model
suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(DBI,warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(stringr,warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(data.table,warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel,warn.conflicts = FALSE)))
###We need to read the peaktable and the data matrices.

###TODO includes dilutionQC.
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
PATH_DB <- args[1]
PATH_DM <- args[2]
TEMP_DM <- args[3]
FOLD_BLANK <- as.numeric(args[4])
QC_FRACTION <- as.numeric(args[5])

dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
all_types <- dbGetQuery(dbb, "SELECT types FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")[, 1]
all_samples <- dbGetQuery(dbb, "SELECT path FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1")[, 1]
dbDisconnect(dbb)



###We get the size of the data matrix.
dm <- fread(PATH_DM,nrows = 2,sep="\t")
cnames <- colnames(dm)
int_prefix <- str_split(cnames,"_")[length(cnames)][[1]][1]
sel_int <- which(str_starts(cnames,int_prefix))
sampled_match_name <- str_match(basename(all_samples),"(.+)\\.[a-zA-Z]+")[,2]
col_match_names <- str_match(cnames[sel_int],paste(int_prefix,"_(.+)\\.csv",sep=""))[,2]

first_line <- 0
batch_size <- ceiling((4000**2)/ncol(dm))

###We get the quantitative 
qcs <- which(tolower(all_types)=="qc")
if(length(qcs)==0) qcs <- which(tolower(all_types)=="sample")
blanks <- tolower(all_types)=="blank"
non_blanks <- seq_along(all_types)
non_blanks <- non_blanks[!blanks]
blanks <- which(blanks)


pqc <- sel_int[qcs]
pblank <- sel_int[blanks]
psample <- sel_int[non_blanks]

num_total <- 0
num_kept <- 0

while(is.data.table(dm)){
  dm <- tryCatch(fread(PATH_DM,skip = first_line,nrows = batch_size),error=function(x){return(NA)},sep="\t")
  num_total <- num_total+nrow(dm)
  if(!is.data.table(dm)) break

  ###Filtering blanks
  if(length(blanks)!=0 & length(non_blanks)!=0){
    int_blank <- apply(dm[,..pblank],1,mean,na.rm=TRUE)
    vfold <- rep(1000,length(int_blank))
    int_sample <- apply(dm[,..psample],1,mean,na.rm=TRUE)
    vfold[is.nan(int_sample)] <- 0
    pok <- (!is.nan(int_blank))|(!is.nan(int_blank))
    vfold[pok] <- int_sample[pok]/int_blank[pok] 
    sel_val <- vfold>=FOLD_BLANK
    dm <- dm[sel_val]
  }
  
  ###Filtering QCs
  if(length(qcs)!=0 & nrow(dm)!=0){
    num_detect <- apply(dm[,..pqc],1,function(x){sum(!is.na(x))})
    frac <- num_detect/length(pqc)
    sel_val <- frac>=QC_FRACTION
    dm <- dm[sel_val]
  }  
  
  num_kept <- num_kept+nrow(dm)
  ###We write the outpu table
  vv <- fwrite(x = dm,file = TEMP_DM,append=TRUE,sep="\t")
  first_line <- first_line+batch_size
}
ff <- file.rename(TEMP_DM, PATH_DM)

