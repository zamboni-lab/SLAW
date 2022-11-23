###Loading the DB connection of samples
# suppressWarnings(suppressMessages(library(RSQLite, warn.conflicts = FALSE)))
# suppressWarnings(suppressMessages(library(DBI, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(stringr, warn.conflicts = FALSE)))
# suppressWarnings(suppressMessages(library(MSnbase, warn.conflicts = FALSE)))
# suppressWarnings(suppressMessages(library(BiocParallel, warn.conflicts = FALSE)))
# suppressWarnings(suppressMessages(library(igraph, warn.conflicts = FALSE)))
# suppressWarnings(suppressMessages(library(rtree, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))

sink(file=stdout())

DEBUG <- F

if (DEBUG==T) {
  #This should be read.
  PATH_DM <- "D:/SW/SLAW_test_data_out/data_full_741d552fefa0759df99c04af0d7f6562.csv"
  PATH_MGF <- "D:/SW/SLAW_test_data_out/spectra_741d552fefa0759df99c04af0d7f6562.mgf"
} else
{
  args <- commandArgs(trailingOnly = TRUE)
  PATH_DM <- args[1]
  PATH_MGF <- args[2]
}

pdm <- fread(PATH_DM,sep="\t")

ISO_DIST_COL <- "isotopic_pattern_abs"
ISO_NAME_COL <- "isotopic_pattern_annot"
MGF_ID_COL <- "slaw_id"
MEAN_INT <- "intensity_mean"
MZ_PRECURSOR <- "mz"
RT_PRECURSOR <- "rt"

MASS_DIFF_REGEXP <- "(?:M\\+([0123456789\\.]+))"
ISO_DIST_REGEXP <- "([0123456789\\.]+)"
MS2_DIST <- "([0-9]+)_\\(e([0123456789\\.]+)\\)"

#We first extract the isotopic informations
mgf_id <- pdm[[MGF_ID_COL]]
fms2 <- which(mgf_id!="")

info_ms2 <- str_match_all(mgf_id[fms2],MS2_DIST)
info_iso_mass <- str_match_all(pdm[[ISO_NAME_COL]][fms2],MASS_DIFF_REGEXP)
info_iso_int <- str_match_all(pdm[[ISO_DIST_COL]][fms2],ISO_DIST_REGEXP)

to_edit <- which((sapply(info_iso_int,nrow)!=0)&(sapply(info_ms2,nrow)!=0))
to_edit_raw_pos <- fms2[to_edit]

mzm <- pdm[fms2,..MZ_PRECURSOR]
rtm <- pdm[fms2,..RT_PRECURSOR]
intm <- pdm[fms2,..MEAN_INT]

##WE build the MS1 spectra for this dataset
iso_spectra <- mapply(info_iso_mass,info_iso_int,as.list(mzm)[[1]],as.list(intm)[[1]],FUN = 
                        function(mziso,intiso,mz,int){
                          if(nrow(mziso)==0){
                            return(matrix(c(mz,int),nrow=1))
                          }
                          intv <- as.numeric(intiso[,2])
                          mzv <- c(mz,mz+as.numeric(mziso[,2]))
                          return(matrix(c(mzv,intv),ncol=2,byrow=FALSE))
             })#iso_spectra


tcon <- file(description = PATH_MGF, open = "a")
writeLines(c('','','# MS1 spectra extracted by SLAW',''),tcon)
## WRITE "MANUALLY" > all packages tested so far disappointed
for (i in 1:nrow(pdm)) {
  spec <- iso_spectra[[i]]
  rows <- c('BEGIN IONS',
            paste0('SCANS=NA'),
            paste0('TITLE=','MS1_slawID=',pdm[i,'slaw_id'],'_mz=',sprintf("%.4f",pdm[i,'mz']), '_RT=',sprintf("%.2f",pdm[i,'rt'])),
            paste0('RTINSECONDS=',pdm[i,'rt']*60),
            paste0('PEPMASS=',pdm[i,'mz']),
            paste0('PRECURSOR_INTENSITY=',as.integer(pdm[i,'intensity_mean'])),
            'MSLEVEL=1',
            paste0('ENERGY=0'),
            paste0('CHARGE=NA'),
            paste0('SLAW_ID=',pdm[i,'slaw_id']),
            paste0('PEAKSCOUNT=',nrow(spec)))
  writeLines(rows,tcon)
  for (j in 1:nrow(spec)) {
    writeLines(paste(sprintf('%.5f',spec[j,1]),spec[j,2]),tcon)
  }
  rows <- c('END IONS','')
  writeLines(rows,tcon)
}
close.connection(tcon)
