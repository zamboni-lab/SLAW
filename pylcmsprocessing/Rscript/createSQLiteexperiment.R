suppressWarnings(suppressMessages(library(optparse,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(DBI,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(RSQLite,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(stringr,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))
suppressWarnings(suppressMessages(library(tools,warn.conflicts = FALSE,quietly = TRUE,verbose = FALSE)))


###CONSTANT
DBSAMPLES <-  "samples" ###Containing peak and path
DBPROCESSINGSAMPLE <-  "processing" ###Containing peak and path
DBPEAKSPICKING <- "peakpicking"
DBNAME <-    "MSexperiment" ###


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
  )

)


opt_parser <-  OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


#####We move to the target experiment repostory.
if (!is.null(opt$directory))
  setwd(opt$directory)

# print(getwd())
####Parsing files and replicates informations
paths <-  NULL
replicates <-  NULL
types <- NULL
if (file.exists(opt$summary)) {
  tsummary <-
    read.table(
      opt$summary,
      header = TRUE,
      stringsAsFactors = FALSE,
      sep = opt$summarysep[[1]]
    )

  ###Getting the pah of the raw files
  if(dir.exists(opt$mzml)){
      paths <- normalizePath(list.files(opt$mzml, full.names = TRUE,pattern=".mzML"))
  }else{
    paths <- normalizePath(list.files(opt$directory, full.names = TRUE,pattern=".mzML"))
  }
  expected_paths <- tsummary[[opt$summarypath]]
  vm <- match(basename(expected_paths), basename(paths))
  if (any(is.na(vm)))
    stop("Missing files :", paste(expected_paths[vm], sep = ", "))
  paths <-  paths[!is.na(vm)]
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
    types <- rep("sample",nrow(tsummary))
  }

} else{
  if(dir.exists(opt$mzml)){
      paths <- normalizePath(list.files(opt$mzml, full.names = TRUE,pattern=".mzML"))
  }else{
    paths <- normalizePath(list.files(opt$directory, full.names = TRUE,pattern=".mzML"))
  }
  replicates <-  rep(1, length(paths))
}

sample_tab <- data.frame(
  id = seq_along(paths),
  path = paths,
  replicate = as.integer(as.factor(replicates)),
  stringsAsFactors = FALSE
)

# message(sample_tab)

#####Writing the sample table and intializing the table


# if (file.exists(opt$database))
#   stop(
#     "Database ",
#     opt$database,
#     " already exists please remove it if you want to repopulat ethe sample table."
#   )

db <- dbConnect(RSQLite::SQLite(), opt$database)
dd <- dbExecute(db,"PRAGMA foreign_keys = ON")

###WE always activate the foreign key constrainct pracma
dd <- dbExecute(db,paste("CREATE TABLE ",DBSAMPLES," (
          id INTEGER PRIMARY KEY,
          path TEXT NOT NULL,
          replicate INTEGER NOT NULL
)",sep=""))


###We then isert the values
dd <- dbAppendTable(conn = db,
             sample_tab,
             name = DBSAMPLES)

####We close the database in every case.
dbDisconnect(db)
