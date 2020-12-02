#' @include common.R


# getHeight <- function(cnames){
#   cnames[which(!is.na(str_match(cnames,"Peak\\.height$")))]
# }
#
# getArea <- function(cnames){
#   cnames[which(!is.na(str_match(cnames,"Peak\\.area$")))]
# }

#' @export
convertPeakTable <- function(peaktable){
  match_names <- sapply(CONVERT_TABLE,function(pat,x){which(str_detect(x,fixed(pat)))},x=colnames(peaktable))
  ###We handle the case of non match thing
  # match_names <- sapply(match_names,function(x){if(is.null(x)){return(NA)};return(x)})
  found_pos <- which(sapply(match_names,length)!=0)
  if(nrow(peaktable)==0){ ##Case of tempy peaktable. Non signal to noise.
    pt <- data.frame()
    for(ic in names(CONVERT_TABLE)[found_pos]){
      pt[[ic]] <- numeric(0)
    }
    return(pt)
  }

  peaktable <-  peaktable[,as.numeric(unlist(match_names[found_pos]))]
  colnames(peaktable) <-  c(names(CONVERT_TABLE)[found_pos])

  ###Supplemtary metric for the diagnosis
  cnames <- colnames(peaktable)
  peakwidth <-  peaktable[,"rt_max"]-peaktable[,"rt_min"]
  assymetry_factor <- (peaktable[,"rt_max"]-peaktable[,"rt"])/(peaktable[,"rt"]-peaktable[,"rt_min"])

  peaktable <- cbind(peaktable,peakwidth,assymetry_factor)
  colnames(peaktable) <- c(cnames,"peakwidth","right_on_left_assymetry")
  return(peaktable)
}
