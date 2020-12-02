setClass("MZmineXML",
         slots = c(
           file = "character",
           list = "list",
           parameters = "list",
           combinations = "data.frame"
         ))

SEPARATOR <- "__"
IGNORED_FIELD <- c("ptype")
P5VAL <- "_EMPTY_"


top5 <- function(val){
  if(is.null(val)||(length(val)==0)){
    return(P5VAL)
  }
  return(val)
}


fromp5 <- function(x=NULL){
  if(x==P5VAL) return(NULL)
  return(x)
}


CONVERT_TABLE <-  list("mz"="row.m.z","rt"="row.retention.time","height"="Peak.height",
                       "intensity"="Peak.area","rt_min"="Peak.RT.start","rt_max"="Peak.RT.end",
                       "mz_min"="Peak.m.z.min","mz_max"="Peak.m.z.max","SN"="Signal.to.Noise")

OUTPUT_MODULES <- data.frame(module=c(
  "__net.sf.mzmine.modules.peaklistmethods.io.csvexport.CSVExportModule__Filename__current_file",
  "__net.sf.mzmine.modules.peaklistmethods.io.gnpsexport.GNPSExportAndSubmitModule__Filename__current_file"
  ),
  name=c("peaktable","msms"),
  format = c("csv","mgf")
,stringsAsFactors = FALSE)
