#' @include MZmineXML.R


############################
###Handling of JSON files###
############################

#' @export
generateTemplateJSON <- function(mxml,file,subset=NULL){
  to_write <-  mxml@parameters
  if(!is.null(subset)) to_write <- to_write[subset]
  to_write <-  toJSON(to_write,pretty=TRUE)
  write(to_write,file)
}

#' @export
importJSON <- function(mxml,file){
  ljson <- fromJSON(file)
  for(nn in names(ljson)){
    mxml@parameters[[nn]] <- ljson[[nn]]
  }
  mxml
}

#############################
###Exporting the Json file###
#############################
