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


path_xml_input <- args[2]

path_xml_output <- normalizePath(args[3])
prefix_xml <- file.path(path_xml_output,"xml")

path_candidates <- args[4]

path_peaktables <- args[5]

path_msms <- args[6]

path_json <- args[7]

message("Json in ",path_json)
if(FALSE){
db <- dbConnect(RSQLite:::SQLite(),dbpath)
}
##The num worker is read form the databse

samplist <- NULL
##TO DEBUBG
if(FALSE){
num_workers <- as.numeric(dbGetQuery(db, "SELECT max_jobs FROM common")[, 1])
}else{
  num_workers <- 1
}


####We get the path form all the sample
samplist <-NULL
if(FALSE){
samplist <- dbGetQuery(db, "SELECT path FROM samples")[, 1]
dbDisconnect(db)
}else{
  samplist <- c("UUU/UUU.mzML","AAA/AAA.mzML")
}


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
generateTemplateJSON(mxml,path_json)
mxml <- importJSON(mxml, path_json)
mxml <- generateParametersCombinationsDataFrame(mxml,maxVal=1)

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

    ###We assign the raw fiel parameter.
    mxml@combinations[, "__net.sf.mzmine.modules.rawdatamethods.rawdataimport.RawDataImportModule__Raw data file names__file"] <-
      rawfile


    PREF_SAMPLE <- paste(prefix,num_rawfile,sep = "_")
    PREF_OUTPUT_PEAKTABLE <- file.path(outp, paste("peaktable_", num_rawfile, sep = ""))
    PREF_OUTPUT_MSMS <- file.path(outp, paste("msms_", num_rawfile, sep = ""))

    vfiles <-
      exportingAllCombinationsOfXML(mxml,
                                    prefix_xml_files = PREF_SAMPLE,
                                    prefix_output = c(PREF_OUTPUT_PEAKTABLE,PREF_OUTPUT_MSMS), verbose = FALSE)
    outfiles <- vfiles[, 2]
    outmsms <- vfiles[,3]

    ###We rename all the files

    ppstemp <-
      data.frame(
        peakpicking = tids,
        sample = rep(num_rawfile, nrow(vfiles)),
        input = vfiles[, 1],
        hash_input = vfiles$hash,
        output_peaktable = outfiles,
        output_msms = outmsms,
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
      tids= 1
    ),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE,
    BPPARAM=bpp
  )

btab <- do.call(rbind, allsampseq)
###We have all the hashes in btab$hash_input
vav <- write.table(btab,file=path_candidates,sep=";",row.names=FALSE)
