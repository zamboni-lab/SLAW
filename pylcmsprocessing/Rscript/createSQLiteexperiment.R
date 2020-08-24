suppressWarnings(suppressMessages(library(optparse,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(DBI,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(stringr,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(tools,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))


###CONSTANT
DBSAMPLES <-  "samples" ###Containing peak and path
DBPROCESSINGSAMPLE <-  "processing" ###Containing peak and path
DBPEAKSPICKING <- "peakpicking"
DBNAME <-    "MSexperiment"


option_list = list(
  make_option(
    c("-d", "--directory"),
    type = "character",
    default = getwd(),
    help = "directory of the experiment",
    metavar = "character"
  ),
  make_option(
    c("-m", "--mzml"),
    type = "character",
    default = "mzML",
    help = "directory of the mzML files [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-s", "--summary"),
    type = "character",
    default = "summary.csv",
    help = "data file summarizing the experiments [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-v", "--summarysep"),
    type = "character",
    default = ",",
    help = "column separation in the summary file [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-r", "--summaryreplicate"),
    type = "character",
    default = "replicate",
    help = "name of the factor column of replicates [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-t", "--summarytype"),
    type = "character",
    default = "type",
    help = "name of the factor column [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-p", "--summarypath"),
    type = "character",
    default = "path",
    help = "name of the summary column of paths [default= %default]",
    metavar = "character"
  ),
  make_option(
    c("-b", "--database"),
    type = "character",
    default = paste("MSexperiment", "sqlite", sep = "."),
    help = "name of default database [default= %default]",
    metavar = "character"
  ),
    make_option(
    c("-o", "--pathms2"),
    type = "character",
    default = NA,
    help = "path to MS2 data if there is some eventually [default= %default]",
    metavar = "character"
  )

)




opt_parser <-  OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
if(FALSE){
  opt <- list(directory="U:/users/Alexis/data/BioMarCoAlaaPos",summary="U:/users/Alexis/data/BioMarCoAlaaPos/summary.csv",
              summarypath="path",summarytype="type",mzml="U:/users/Alexis/data/BioMarCoAlaaPos")
  
}


#####We move to the target experiment repostory.
if (!is.null(opt$directory))
  setwd(opt$directory)

# print(getwd())
####Parsing files and replicates informations
paths <-  NULL
replicates <-  NULL
types <- NULL
MSlevels <- NULL

if (file.exists(opt$summary)) {
  tsummary <-
    read.table(
      opt$summary,
      header = TRUE,
      stringsAsFactors = FALSE,
      sep = opt$summarysep[[1]]
    )

  path_dir <- NULL
  ###Getting the pah of the raw files
  if(dir.exists(opt$mzml)){
    path_dir <- opt$mzml
  }else{
    path_dir <- opt$directory
  }
  paths <- normalizePath(list.files(path_dir, full.names = TRUE,pattern=".mzX?ML"))

  expected_paths <- tsummary[[opt$summarypath]]
  vm <- match(basename(expected_paths), basename(paths))
  if (any(is.na(vm)))
    stop("Missing files :", paste(expected_paths[is.na(vm)], sep = ", "))
  ##We reorder path
  
  if (anyDuplicated(vm)){
    stop("Duplicated paths in summary.txt:",paste(expected_paths[duplicated(vm)], sep = ", "))
  }
  
  
  paths <- paths[vm]
  rnames <- opt$summary
  if (opt$summaryreplicate %in% colnames(tsummary)) {
    replicates <-  tsummary[[opt$summaryreplicate]]
  } else{
    replicates <-  rep(1, length(paths))
  }
  
  if (opt$summarytype %in% colnames(tsummary)){
    types <- tsummary[[opt$summarytype]]


    ###Three type of terminology are authorized, QC,blank,sample and QC and a number.
    msamp <- types=="sample"
    mQC <- types=="QC"
    mBlank <- types=="blank"

    ###We extrct the dilution QCs.
    vmatch <-str_match(types,"QC([0-9]+)")
    mDilQCs <- !is.na(vmatch)

  }else{
    types <- rep("sample",length(paths))
  }

} else{
  if(dir.exists(opt$mzml)){
      paths <- normalizePath(list.files(opt$mzml, full.names = TRUE,pattern=".mzX?ML"))
  }else{
    paths <- normalizePath(list.files(opt$directory, full.names = TRUE,pattern=".mzX?ML"))
  }
  types <- rep("sample",length(paths))
  replicates <-  rep(1, length(paths))
}

levels_str <- rep("MS2",length(paths))
ms1_str <- rep("MS1",length(paths))
levels_str <- ifelse(types=="MS2",levels_str,ms1_str)

sample_tab <- data.frame(
  id = seq_along(paths),
  path = paths,
  level = levels_str,
  types = types,
  replicate = as.integer(as.factor(replicates)),
  stringsAsFactors = FALSE
)

###We check if an MS2 argument is present if it the case we add the samples in the table as MS2
if(!is.na(opt$pathms2)){
  ##We add the MS2
  path_ms2 <- list.files(opt$pathms2,full.names=TRUE)
  ##Now we rebuild the next step
  supp_tab <- data.frame(
    id = seq_along(path_ms2)+nrow(sample_tab),
    path = path_ms2,
    level = rep("MS2",length(path_ms2)),
    types = rep("MS2",length(path_ms2)),
    replicate = as.integer(rep(0,length(path_ms2))),
    stringsAsFactors = FALSE
  )
  sample_tab <- rbind(sample_tab,supp_tab)
}

db <- dbConnect(RSQLite::SQLite(), opt$database)
dd <- dbExecute(db,"PRAGMA foreign_keys = ON")

###WE always activate the foreign key constrainct pracma
dd <- dbExecute(db,paste("CREATE TABLE IF NOT EXISTS ",DBSAMPLES," (
          id INTEGER PRIMARY KEY,
          path TEXT NOT NULL,
          level TEXT NOT NULL,
          types TEXT NOT NULL,
          replicate INTEGER NOT NULL
)",sep=""))


###We then isert the values
dd <- dbAppendTable(conn = db,
             sample_tab,
             name = DBSAMPLES)

####We close the database in every case.
dbDisconnect(db)
