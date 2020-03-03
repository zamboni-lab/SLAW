###A set of model peaks

###Model has a relative value

setClass("LCMSAlignerReferences",slots=list(peaks="data.frame",
                                            idxClasses="list",
                                            parameters="list",bpp = "ANY"))

setClass("LCMSAlignerFileHandler",
         slots = list(files="environment",inv_files="environment",
                      id_table="environment",num_files="numeric"))

setClass("LCMSAlignerDataStorage",
         slots = list(dir="character",files="list",frequency="numeric"))


#' LCMSAlignerModel
#'
#' @slot references The mains paramters which are reused. 
#' @slot peaks A dataframe storing informations about ALL the peaks at any stage.
#' @slot data A list of all storing all the occurences of peaks currently in memoery
#' @slot files LCMSAlignerFileHandler The file handler system ensuring that multiple files are not processed multiple files
#' @slot storage LCMSAlignerDataStorage A file system storing all the data 
#' @slot rt_correction A list storing all the corrected RTs just for fun.
#' @slot files_in_memory A list storing the ids of all the character in memory.
#'
#' @return
#' @export
#'
#' @examples
#' print("Examples to be put here")
setClass("LCMSAlignerModel",slots=list(references="LCMSAlignerReferences",
                                       peaks="data.frame",
                                       data="matrix",
                                       files="LCMSAlignerFileHandler",
                                       storage = "LCMSAlignerDataStorage",
                                       rt_correction="list",
                                       files_in_memory = "numeric",
                                       save_interval = "numeric",
                                       order_peaks = "numeric",
                                       clustering = "numeric",
                                       max_id = "numeric"))
