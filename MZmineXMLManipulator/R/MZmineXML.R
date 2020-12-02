#' @include common.R
#' @include XMLhandler.R


#' @export
MZmineXML <- function(path) {
  mxml <- new("MZmineXML")
  mxml@file <-  path
  mxml@list <-  parseMzMineXML(path)
  leav <- leaves(mxml, usedOnly = FALSE)
  mxml@parameters <-  sapply(leav,getLeaf,mxml=mxml,simplify = FALSE)
  # mxml@parameters <- vector(mode="list",length = length(leav))
  names(  mxml@parameters) <- leav
  mxml
}


#' @export
setMethod("show", "MZmineXML", function(object) {
  cat("An MZmineXML object generated from", object@file,"containing",length(object@parameters),"parameters.\n")
  if(nrow(object@combinations)!=0)
  cat("It contains",nrow(object@combinations),"rows.\n")
})

#' @export
getLeaf <- function(mxml, path) {
  path <-  str_split(path, pattern = fixed(SEPARATOR))[[1]][-1]
  recur_access(mxml@list, path)
}

getParam <- function(mxml,param){
  mxml@parameters[[param]]
}


recur_access <- function(terms, path) {
  if (length(path) == 0) {
    return(terms)
  }

  if(length(path)==2 && !suppressWarnings(is.na(as.numeric(path[2])))){
    return(terms[[as.numeric(path[2])]])
  }
  return(recur_access(terms[[path[1]]], path[-1]))
}

########################
###Parameters setting###
########################
#' @export
assignParameter <- function(mxml,path,value){
  value <-  top5(value)
  posPar <- match(path,names(mxml@parameters))
  if(is.na(posPar)) stop("Unknown parameter ",path," \n")
  if(all(getParam(mxml,path)==value)) return(mxml)
  ssplit <-  str_split(path,pattern = fixed(SEPARATOR))[[1]][-1]
  mxml@list <- recur_assign(mxml@list,path=ssplit,value=value)
  mxml@parameters[[path]] <- value
  mxml
}

recur_assign <- function(terms, path, value) {
  if (length(path) == 0) {
    return(value)
  }

  if(length(path)==2 && !suppressWarnings(is.na(as.numeric(path[2])))){
    terms[[as.numeric(path[2])]] <- value
    return(terms)
  }

  terms[[path[1]]] <- recur_assign(terms[[path[1]]], path[-1], value)
  return(terms)
}








################################
###Arguments extractions path###
################################

filterUseless <- function(x) {
  return(x[!(names(x) %in% IGNORED_FIELD)])
}

isLeaf <- function(x) {
  !any(sapply(x, is.list))
}


isUsed<- function(x) {
  isUsedLeaf(x)
  # all(!sapply(x,isUsedLeaf))
}

isUsedLeaf <- function(x) {
  if (length(x) == 0) {
    return(FALSE)
  } else{
    return(TRUE)
  }
}

isModule <- function(x){
  if( !( "ptype" %in% names(x))) return(FALSE)
  x$ptype==4
}

isParam6 <- function(x){
  if( !( "ptype" %in% names(x))) return(FALSE)
  x$ptype==6
}

treat_duplicates <-  function(x){
  counter <- 1
  in_seq <-  FALSE
  for(i in seq_along(x)){
    if(i==1) next
    if(x[i-1]==x[i]){
      x[i-1] <- paste(x[i-1],counter,sep=SEPARATOR)
      counter <- counter+1
      in_seq <-  TRUE
    }else{
      if(in_seq){
        in_seq <-  FALSE
        x[i-1] <- paste(x[i-1],counter,sep=SEPARATOR)
      }
      counter <- 1
    }
  }
  return(x)
}


leaves <-  function(mxml, usedOnly = FALSE) {
  res <- dfs_trav_prefix(mxml@list, "",usedOnly,depth=0)
  res <- unname(unlist(res, recursive = TRUE))
  res <- res[!is.na(res)]
  temp <- treat_duplicates(res)
  temp[!is.na(temp)]
}


dfs_trav_prefix <- function(node, prefix,depth, usedOnly = FALSE){

  if(typeof(node)=="character"){
    return(prefix)
  }

  ####special case if it s a module parameter we only rocess the correct module
  if((depth>=2)&&(isModule(node))){
    psel <- match(node$selected,names(node))
    return(dfs_trav_prefix(node[[psel]], paste(prefix, names(node)[psel],sep=SEPARATOR), depth+1, usedOnly = usedOnly))
  }

  ###We check if the type of parameters, if they are not selected we don t consider them.
  if((depth>=2)&&(isParam6(node))){
    if(node$selected){
      node <- filterUseless(node)
      ##We parse all the separator
      return(mapply(node, paste(prefix, names(node), sep = SEPARATOR), FUN =
               dfs_trav_prefix,depth=depth+1,usedOnly=usedOnly))
    }else{
      return(NA_character_)
    }
  }

  node <- filterUseless(node)
  if (isLeaf(node)) {
    if(usedOnly && !isUsed(node)){
      return(NA_character_)
    }else{
      temp <- filterUseless(node)
      return(c(paste(prefix, names(temp), sep = SEPARATOR)))
    }
  }
  return(mapply(node, paste(prefix, names(node), sep = SEPARATOR), FUN =
                  dfs_trav_prefix,depth=depth+1,usedOnly=usedOnly))
}

##############################################
###Generating all combination of parameters###
##############################################

#' @export
generateParametersCombinationsDataFrame <- function(mxml,maxVal=200,seed=512){
  mxml@combinations <- do.call("expand.grid",c(mxml@parameters,stringsAsFactors=FALSE))
  if(!is.null(maxVal) & nrow(mxml@combinations)>maxVal){
    set.seed(seed)
    mxml@combinations <- mxml@combinations[sample(1:nrow(mxml@combinations),size = maxVal),,drop=FALSE]
  }
  mxml
}

#' @export
setParametersCombinationsDataFrame <- function(mxml,df){
  mxml@combinations <- df
  return(mxml)
}



#' @export
getCombinations <-  function(mxml){
  mxml@combinations
}

#' @export
setCombinations <-  function(mxml,value){
  mxml@combinations <-  value
  mxml
}


###########################################
###Exporting all the combinations of XML###
###########################################

#' @export
writeXML <-  function(lmxml,file){
  output_xml(lmxml, file = file, compression = 0, indent = TRUE,
             indentchar="\t", prefix = "<?xml version=\"1.0\"?>\n",
             doctype = NULL, encoding = getEncoding(doc))
}



getOutputPrefix <-  function(mxml,prefix){
  pout <- match(names(mxml@parameters),OUTPUT_MODULES[,1])
  pok <-  which(!is.na(pout))

    ###Now we give the
  if(length(pok)==0){
    return(NULL)

  }else if(is.na(prefix)){
    stop("Output filed module spotted:",names(mxml@parameters)[pout[pok]]," you need to give an output prefix.")
  }
  return(pout)
}



hash_function <- function(x){
  return(md5(x))
}


#' @export
exportingAllCombinationsOfXML <-  function(mxml,prefix_xml_files=NULL,summary_file=NULL,prefix_output=NA_character_,
                                           hash_checking=TRUE,verbose=TRUE){
  if(is.null(prefix_xml_files)&is.null(summary_file)) stop("No file to process")
  if(verbose)
  cat("Generating XML batch files for MZmine.")
  fnames <- paste(prefix_xml_files,paste(1:nrow(mxml@combinations),".xml",sep=""),sep="_")

  ### We change the names of the output to reflect the names multiple combinations of mxml
  vout <- getOutputPrefix(mxml,prefix_output[1])
  pok <-  NULL
  if(!is.null(vout)){
    pok <-  vout[which(!is.na(vout))]
    if(length(prefix_output)!=length(pok)){
      stop("EE: ",paste(prefix_output,collapse=",")," EE: ",paste(pok,collapse=","))
#multiple output found, please specify another output directory.")
      stop("multiple output found, please specify another output directory.")
    }
    l_out <- list()
    for(i in pok){
      l_out[[OUTPUT_MODULES[i,"name"]]] <- character(nrow(mxml@combinations))
    }
  }

  names_supp_summary <- paste("output",OUTPUT_MODULES[pok,"name"],sep="_")

  hash_vals <- NULL
  if(hash_checking){
    hash_vals <- character(length(fnames))
  }
  pb <- NULL
  if(verbose)
  pb <- txtProgressBar(style=3)
  if(!is.null(prefix_xml_files)){
    for(i in 1:nrow(mxml@combinations)){
      tmxml <- mxml
      if(verbose)
      setTxtProgressBar(pb, i/nrow(mxml@combinations))
      for(parv in colnames(tmxml@combinations)){
        tmxml <- assignParameter(mxml = tmxml,path = parv,value = mxml@combinations[i,parv])
      }
      hash_val <- NULL
      f_exists <- FALSE
      if(hash_checking){
        hash_val <-hash_function(paste(c(as.character(mxml@combinations[i,]),as.character(i)),collapse = "_"))
        hash_vals[i] <- hash_val
        fnames[i] <- paste(prefix_xml_files,"_",hash_val,".xml",sep="")
      }

      ####We hash the parameter and always include the name of the file
      if(!is.null(vout)){
        for(idj in seq_along(pok)){
          j <- pok[idj]
          output_name <- paste(prefix_output[idj],"_",hash_val,".",
                             OUTPUT_MODULES[j,"format"],sep="")
          if(file.exists(output_name) & hash_checking) f_exists <- TRUE
          tmxml <- assignParameter(mxml = tmxml,path = OUTPUT_MODULES[j,"module"],value = output_name)
          l_out[[OUTPUT_MODULES[j,"name"]]][i] <- output_name
        }
      }

      if(!f_exists){
        temp <- createMzMineXML(tmxml@list)
        output_xml(temp,fnames[i])
      }
    }
  }
    onames <- names(mxml@combinations)
    summary_df <-  NULL

    largs <- list(fnames)
    lnames <- c("file")

    if(!is.null(vout)){
      df_temp <- data.frame(out1=l_out[[1]],stringsAsFactors = FALSE)
      if(length(l_out)>1){
        for(j in 2:length(l_out)){
          df_temp[[paste("output",j,sep="")]] <- l_out[[j]]
        }
      }
      largs <- append(largs,df_temp)
      lnames <- append(lnames,names_supp_summary)
    }

    if(hash_checking){
      largs <- append(largs,data.frame(hash=hash_vals,stringsAsFactors = FALSE))
      lnames <- append(lnames,"hash")
    }

    largs <- append(largs,mxml@combinations)
    lnames <- append(lnames,onames)

    summary_df <- as.data.frame(do.call(cbind,largs))

    colnames(summary_df) <- lnames

  if(!is.null(summary_file)){
    write.table(file = summary_file,x = summary_df,sep=";",row.names = FALSE)

  }
    return(invisible(summary_df))
}


