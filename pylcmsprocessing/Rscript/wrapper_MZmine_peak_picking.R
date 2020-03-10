suppressMessages(library(DBI))
suppressMessages(library(RSQLite))
suppressMessages(library(BiocParallel))
suppressMessages(library(jsonlite))
suppressMessages(library(MZmineXMLManipulator))
suppressMessages(library(tools))
suppressMessages(library(stringr))


args <- commandArgs(trailingOnly = TRUE)

###We can get everything from the parameter
dbpath <-  args[1]

path_json <- args[2]

path_xml_input <- args[3]

path_xml_output <- normalizePath(args[4])
prefix_xml <- file.path(path_xml_output,"xml")

path_summary <- args[5]

path_summary_templates <- args[6]

path_templates <- normalizePath(args[7])
prefix_templates <- file.path(path_templates,"templates")

path_candidates <- args[8]

path_peaktables <- normalizePath(args[9])

peak_picking_id <- as.numeric(args[10])
db <- dbConnect(RSQLite:::SQLite(),dbpath)
##The num worker is read form the databse

num_workers <- as.numeric(dbGetQuery(db, "SELECT max_jobs FROM common")[, 1])

####We get the path form all the sample
samplist <- dbGetQuery(db, "SELECT path FROM samples WHERE level='MS1'")[, 1]
message(paste(samplist,sep="|"))

dbDisconnect(db)

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

bpp <- NULL
if(get_os()=="win"){
  bpp <- SnowParam(workers = num_workers)
}else{
  bpp <- MulticoreParam(workers = num_workers)
}

####SELECT AT MOST MAX PARAMETERS, THE SAME FOR EVERY FILE, GENERATE THE ASSOCIATED XML PUT THEM INTO DB
mxml <- MZmineXML(path_xml_input)

###We first sample the parameter space
pjson <-  file_path_as_absolute(path_json)
mxml <- importJSON(mxml, pjson)


mxml <- generateParametersCombinationsDataFrame(mxml,maxVal=2)

summary_table <- read.table(path_summary,sep=";",stringsAsFactors=FALSE,header=TRUE,comment.char="",check.names=FALSE)


mxml <-  setCombinations(mxml,summary_table)

####We write all the selected parameter in a directory

templates <-
  exportingAllCombinationsOfXML(
    mxml,
    summary_file = path_summary_templates,
    prefix_xml_files = prefix_templates,
    prefix_output = "NOTUSED"
  )
templates <- templates[, 1]

####We read the name of the file
####Now we write the sleected cobination in a parameter
singleXMLgeneration <-
  function(rawfile,
           num_rawfile,
           prefix,
           mxml,
           outp, tids) {
    suppressMessages(library(MZmineXMLManipulator))
    suppressMessages(library(tools))

    getRawName <- function(x){
      return(str_split(str_sub(basename(x)),fixed("."))[[1]][1])
    }

    ###We assign the raw fiel parameter.
    mxml@combinations[, "__net.sf.mzmine.modules.rawdatamethods.rawdataimport.RawDataImportModule__Raw data file names__file"] <-
      rawfile


    PREF_SAMPLE <- paste(prefix,num_rawfile,sep = "_")
    PREF_OUTPUT <- file.path(outp, paste("peaktable_", num_rawfile,"_",getRawName(rawfile), sep = ""))

    vfiles <-
      exportingAllCombinationsOfXML(mxml,
                                    prefix_xml_files = PREF_SAMPLE,
                                    prefix_output = PREF_OUTPUT, verbose = FALSE)
    outfiles <- vfiles[, 2]

    ###We rename all the files

    ppstemp <-
      data.frame(
        peakpicking = tids,
        sample = rep(num_rawfile, nrow(vfiles)),
        input = vfiles[, 1],
        hash_input = vfiles[,3],
        output = outfiles,
        stringsAsFactors = FALSE
      )
    return(ppstemp)
  }

allsampseq <-
  bpmapply(
    FUN = singleXMLgeneration,
    samplist,
    seq_along(samplist),
    MoreArgs = list(
      mxml = mxml,
      prefix = prefix_xml,
      outp = path_peaktables,
      tids= 1:nrow(summary_table)
    ),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE,
    BPPARAM=bpp
  )

btab <- do.call(rbind, allsampseq)
###We have all the hashes in btab$hash_input
vav <- write.table(btab,file=path_candidates,sep=";",row.names=FALSE)
