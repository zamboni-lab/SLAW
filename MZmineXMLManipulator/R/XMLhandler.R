#' @include common.R

###Extracting a json from the mzmine parameter>


#
# Multiple cases :
#   If it is a batch step dictionnary all children in dictionnary
#
# FIRST LEVELS ALWAYS REPRESENT A BATCH STEP:
#   <batchstep method="net.sf.mzmine.modules.rawdatamethods.peakpicking.massdetection.MassDetectionModule">
#
#   SECOND LEVEL IS AWAYS A SET OF PARAMETERS:
#   1st type : <parameter name="Raw data files" type="BATCH_LAST_FILES"/>
#   2nd type : <parameter name="Mass list">masses</parameter>
#   3nd type : <parameter name="m/z tolerance">
#   <absolutetolerance>0.003</absolutetolerance>
#   <ppmtolerance>15.0</ppmtolerance>
#   </parameter>
#   4 th type : <parameter name="Mass detector" selected="Centroid">
#   In this case the parameter is followd by a list of module.
#
# Thid elvel if there is any module parameters.
# <module name="Centroid">
#   <parameter name="Noise level">100.0</parameter>
#   </module>
#   <parameter name="Noise level"/>
# 5th type <parameter name="Min m/z peak width"/>


########################
####FROM XML TO List####
########################

parseMzMineXML <-  function(xmlf){
  rtree <- xmlToList(xmlParse(xmlf))
  tempres <- lapply(rtree, parseBatchStep)
  temp_names <- sapply(tempres,'[[',i="name")

  ###We check for duplicated batch step
  if(any(duplicated(temp_names))){
    temp_names[duplicated(temp_names)] <- paste(temp_names[duplicated(temp_names)],
                                                1:sum(duplicated(temp_names)),sep="")
  }

  names(tempres) <- temp_names
  tempres <- sapply(tempres,'[[',i="val")
  return(tempres)
}


parameterType <-  function(pa){
  if(length(pa)==1 && names(pa)=="name") return(5)
  if(typeof(pa)=="list"){
    if((length(pa)==2)&&(all(names(pa)==c("text",".attrs")))){
      return(2)
    }
    if("module" %in% names(pa)) return(4)
    if("selected" %in% names(pa$.attrs)) return(6)
    if(length(pa)>=2){
      return(3)
    }
  }
  if(typeof(pa)=="character"){
    return(1)
  }
}


parseBatchStep <-  function(bs){
  sizeMax <-  length(bs)
  ###We determine the type of the parameter :
  ptype <-sapply(bs[1:(sizeMax-1)],parameterType)
  vpeaks <- mapply(bs[1:(sizeMax-1)],ptype,FUN=parseParameters,SIMPLIFY = FALSE)
  vnames <- sapply(vpeaks,'[[',i="name")
  rlist <- sapply(vpeaks,'[[',i="val",simplify = FALSE)
  names(rlist) <-  vnames
  name <- bs$.attrs[[1]]
  return(list(name=name,val=rlist))
}


parseParameters <-  function(pa,type){
  temp <- NULL
  if(type==1){
    temp <- parseParameters_1(pa)
  }else if(type==2){
    temp <- parseParameters_2(pa)
  }else if(type==3){
    temp <- parseParameters_3(pa)
  }else if(type==4){
    temp <- parseParameters_4(pa)
  }else if(type==5){
    temp <- parseParameters_5(pa)
  }else if(type==6){
    temp <- parseParameters_6(pa)
  }
  return(temp)
}


coerceParameter <- function(x){
  if(length(x)==0) return(NULL)
  if(!suppressWarnings(is.na(as.numeric(x)))) return(as.numeric(x))
  return(x)
}


parseParameters_1 <-  function(x){
  ###We elways try to extract hte name as the place with the name of 1
  pname <- match("name",names(x))
  if(is.na(pname)){
    name <- x[[1]]
    pname <- 1
  }else{
    name <- x[[pname]]
  }

  if(length(x)==1) return(list(name=x[[1]],val=list(ptype=1)))
  x <- x[-pname]
  vnames <- names(x)
  vpar <-  sapply(x,coerceParameter,simplify = FALSE)
  names(vpar) <- vnames
  list(name=name,val=c(vpar,ptype=1))
}

parseParameters_2 <-  function(x){
  templist <- as.list(x$.attrs)
  name <-  templist$name
  templist$name <- NULL
  templist$text <- coerceParameter(x$text)
  list(val=c(templist,ptype=2),name=name)
}

parseParameters_3 <-  function(x){
  num_val <- length(x)
  vnames <- names(x)[1:(num_val-1)]
  name <- x$.attrs
  res_list <- x[1:(num_val-1)]
  res_list <-  sapply(res_list, coerceParameter,simplify = FALSE)
  res_list <- as.list(res_list)
  names(res_list) <-  vnames
  list(val=c(res_list,ptype=3),name=name)
}

parseModule <-  function(x){
  if(typeof(x)=="character"){
    return(list(name=x["name"],val=top5(x[-1])))
  }
  ###If it a one line module without arguments
  if( (length(x)==1) && (names(x)=="name")){
    return(list(name=x[[1]],val=list()))
  }

  ###A module is similar ot a bathc step
  tempres <- parseBatchStep(x)
  return(tempres)
}

parseParameters_4 <-  function(x){
  num_val <- length(x)
  res_modules <- lapply(x[1:(num_val-1)],parseModule)
  vnames <- sapply(res_modules,'[[',i="name")
  rlist <- lapply(res_modules,'[[',i="val")
  names(rlist) <-  vnames
  sel <- unname(x$.attrs["selected"])
  return(list(val=c(rlist,selected=sel,ptype=4),name=x$.attrs["name"]))
}


parseParameters_5 <-  function(x){
  temp <- list(val=list(ptype=5),name=x["name"])
  return(temp)
}


###This is an optional arg
parseParameters_6 <- function(x){
  ##The last argument alway includes the name and everything else.
  ##The rest should be parsed as standard argument
  sizeMax <- length(x)
  ptype <-sapply(x[1:(sizeMax-1)],parameterType)
  res_params <- mapply(x[1:(sizeMax-1)],ptype,FUN=parseParameters,SIMPLIFY = FALSE)
  vnames <- sapply(res_params,'[[',i="name")
  rlist <- lapply(res_params,'[[',i="val")
  names(rlist) <-  vnames
  sel <- unname(x$.attrs["selected"])
  return(list(val=c(rlist,selected=sel,ptype=6),name=x$.attrs["name"]))
}

###################
####List to XML####
###################

#
# newXMLNode("block", "xyz", attrs = c(id = "bob"),
#            namespace = "fo",
#            namespaceDefinitions = c("fo" = "http://www.fo.org"))
#
#

createMzMineXML <- function(dl){

  root <- newXMLNode("batch")

  ###We remove all the number of the names
  vnames <- str_extract(names(dl),pattern = "[a-zA-Z\\.]+")


  ####We need to process all the things at the batch step levels
  bs <- mapply(dl,vnames,FUN=createBatchStepXMLNode,SIMPLIFY = FALSE)

  addChildren(root,as.list(bs))
  return(root)
}

createBatchStepXMLNode <- function(bs,name){

  params <- mapply(bs,names(bs),FUN=createXMLParamNode,SIMPLIFY = FALSE)
  newXMLNode("batchstep",attrs = c(method=name),.children = params)
}

createXMLParamNode <- function(param,name){
  ###We extact the type of the parameters.
  type <- param$ptype
  param <- param[1:(length(param)-1)]
  node <- createParameterXMLnodeType(param,name,type)

}


createParameterXMLnodeType <- function(param,name,type){
  if(type==1){
    return(createXMLParamNode1(param,name))
  }
  if(type==2){
    return(createXMLParamNode2(param,name))
  }
  if(type==3){
    return(createXMLParamNode3(param,name))
  }
  if(type==4){
    return(createXMLParamNode4(param,name))
  }
  if(type==5){
    return(createXMLParamNode5(param,name))
  }
  if(type==6){
    return(createXMLParamNode6(param,name))
  }
}


createXMLParamNode1 <- function(x,name){
  #<parameter name="Peak lists" type="BATCH_LAST_PEAKLISTS"/>
  return(newXMLNode("parameter",attrs = c(x,"name"=name)))
}


createXMLParamNode3 <- function(x,name){
  # <parameter name="Peak duration range (min)">
  #   <min>0.0</min>
  #   <max>10.0</max>
  #   </parameter>
  alln <- mapply(x,names(x),FUN=function(x,y){
    if(length(x)==0) return(newXMLNode(y))
    return(newXMLNode(y,unname(x)))
  })
  newXMLNode("parameter",attrs = list(name=name),.children = alln)
}

createXMLParamNode2 <- function(x,name){
  # <parameter name="Mass list name">masses</parameter>
  textv <- x$text
  supp_attr <- NULL
  if(length(x)>1){
    supp_attr <- names(x)
    psel <- supp_attr!="text"
    supp_attr <- x[psel]
    names(supp_attr) <- names(x)[psel]
    supp_attr <- c("name"=name,supp_attr)
  }else{
    supp_attr <- c("name"=name)
  }
  return(newXMLNode("parameter",textv,attrs = supp_attr))
}

createXMLParamNode4 <- function(x,name){
  # <parameter name="Mass detector" selected="Centroid">
  sel <- x[[length(x)]]
  x <- x[1:(length(x)-1)]
  mods <-  mapply(x,names(x),FUN=createXMLModule)
  newXMLNode("parameter",attrs = c("name"=name,"selected"=sel),.children = mods)

}

createXMLParamNode5 <- function(x,name){
  # <parameter name="Noise level"/>
  return(newXMLNode("parameter",attrs = c("name"=name)))
}

createXMLParamNode6 <- function(x,name){
  # <parameter name="Mass detector" selected="Centroid">
  sel <- x[[length(x)]]
  x <- x[1:(length(x)-1)]
  mods <-  mapply(x,names(x),FUN=createXMLParamNode)
  newXMLNode("parameter",attrs = c("name"=name,"selected"=sel),.children = mods)

}



createXMLModule <- function(x,name){
  # <module name="Local maxima">
  if(length(x)!=length(names(x))) return(newXMLNode("module",attrs = c("name"=name)))
  pars <-  mapply(x,names(x),FUN=createXMLParamNode)
  newXMLNode("module",attrs = c("name"=name),.children = pars)
}

####Internal method is always called


output_xml <-function (doc, file = NULL, compression = 0, indent = TRUE,
                       indentchar="\t", prefix = "<?xml version=\"1.0\"?>\n",
                       doctype = NULL, encoding = NULL, ...){
  # browser()
  ####We first erase the old file
  if(file.exists(file)) file.remove(file)
  if(!is.null(encoding)){
    encoding <-  XML::getEncoding(doc)
  }


  f <- file(file,open="a")
  if (is.na(encoding) || length(encoding) == 0 || encoding ==
      "")
    encoding = character()
  ans = .Call("RS_XML_printXMLNode", doc, as.integer(0), as.integer(indent),
              as.logical(indent), as.character(encoding), XML:::getEncodingREnum(as.character(encoding)),
              PACKAGE = "XML")
  if(nchar(prefix)>0){
    cat(prefix, file = f)
  }
  if (length(file)) {
    cat(ans, file = f)
    close.connection(f)
  }
  else ans
}
