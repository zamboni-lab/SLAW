###Given as input a list of file or a list of ids, extract the associated EICs.
library(DBI)
library(RSQLite)


ars <- commandArgs(trailingOnly = TRUE)
args <- c(
  "U:/users/Alexis/sandbox/input_features.txt",
  "U:/Alexis/sandbox/input_rawfiles.txt",
  "U:/users/Alexis/examples_lcms_workflow/output/processing_db.sqlite",
  "",
  "",
  "0.05",
  "0.02",
  "3"
)


INPUT_FEATURES <- args[1]  
INPUT_RAW_FILES <- args[2]
PATH_DB <- args[3]
OUTPUT_HDF5 <- args[4]
OUTPUT_PDF <- args[5]
TOL_RT <- as.numeric(args[6])
TOL_MZ <- as.numeric(args[7])
NCORES <- args[8]

MARGIN_RT <- 3/60 #(3s)
MARGIN_MZ <- 0.007#MZ tolerance


####We check the class of the input
input_signals <- read.table(INPUT_FEATURES,sep=",",header=FALSE,stringsAsFactors = FALSE)


###We read the assocaited data.frame
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
PATH_PEAKTABLE <- dbGetQuery(dbb, "SELECT peaktable FROM PEAKPICKING")[1, 1]
dbDisconnect(dbb)


dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
raw_files <- dbGetQuery(dbb, "SELECT path FROM samples")[, 1]
dbDisconnect(dbb)


if(FALSE){
  library(stringr)
  PATH_PEAKTABLE <- str_replace(string = PATH_PEAKTABLE,pattern = fixed("/sauer1"),replacement = "U:")
  raw_files <- str_replace(string = raw_files,pattern = fixed("/sauer1"),replacement = "U:")
}

vmz <- NULL
vrt <- NULL

peaktable <- read.table(PATH_PEAKTABLE,sep=",",header=TRUE)

if(ncol(input_signals)==1){
  ###Then we consider that the data have been directly removed.
  sel_idx <- input_signals[,1]
  vmz <- peaktable[input_signals[,1],1]
  vrt <- peaktable[input_signals[,1],2]
}else if(ncol(input_signals)>=2){
    ###authorized tilte are m rt and title.
    ###If this case we need to match the features
    sel_idx <- matchFeatures(input_signals[,1],inputs_signals[,2])
    
    vmz <- input_signals
    vrt <- input_signals[,2]
    if(any(vrt>75)){
      ###we convert the retention time to minutes
      vrt <- vrt/60
    }

}



title <- rep("",length(vmz))
if(ncol(input_signals)>=3){
  if("title" %in% colnames(input_signals)){
    title <- input_signals[,"title"]
  }else{
    title <- input_signals[,3]
  }
}

supp_title <- paste("RT: ",sprintf("%0.2f",vrt-MARGIN_RT),"-", sprintf("%0.2f",vrt+MARGIN_RT),
               "\nMZ:",sprintf("%0.4f",vmz-MARGIN_MZ),"-",sprintf("%0.4f",vmz+MARGIN_MZ),sep="")
title <- paste(title,supp_title,sep=" ")

####
extract_EIC_raw_file <- function(path_xraw,path_peaktable,mz,rt,mz_margin,rt_margin){
  ###We also map the peaktbale in the data 
  library(xcms)
  library(igraph)
  
  if(!file.exists(path_peaktable)) stop(paste(path_peaktable,"file does not exist"))
  peaktable <- read.table(path_peaktable,header=TRUE,sep=",")
  if(!file.exists(path_xraw)) stop(paste(path_xraw,"file does not exist"))
  xraw <- xcmsRaw(path_xraw)
  
  matchLCMSsignals <-
    function(mz_data,
             rt_data,
             mz_ref,
             rt_ref,
             tol_mz=0.01,
             tol_rt = 10){
      ###We map the function to the napping
      mode <-  match.arg(mode)
      # vm <- mineMS2:::matchMzs(mz_ref, mz_data, ppm = ppm,dmz=dmz)
      vm <- xcms:::fastMatch(mz_ref,mz_data,tol = tol_mz)
      pf <-  which(!sapply(vm, is.null))
      
      df_edges <-
        mapply(
          vm[pf],
          pf,
          seq_along(pf),
          FUN = function(x,
                         xi,
                         yi,
                         mz_data,
                         rt_data,
                         mz_ref,
                         rt_ref,
                         alpha_rt,
                         intlog) {
            if (is.na(x)) {
              return(data.frame(
                from = numeric(0),
                to = numeric(0),
                weight = numeric(0),
                mz1=numeric(0),
                rt1=numeric(0),
                mz2=numeric(0),
                rt2=numeric(0)
              ))
            } else{
              ###deviation in ppm
              cmz <- abs(mz_ref[xi] - mz_data[x])/ tol_mz
              crt <- abs(rt_ref[xi] - rt_data[x]) / tol_rt
              if(all(crt>(2*tol_rt)))
                return(data.frame(
                  from = numeric(0),
                  to = numeric(0),
                  weight = numeric(0),
                  mz1=numeric(0),
                  rt1=numeric(0),
                  mz2=numeric(0),
                  rt2=numeric(0)
                ))
              
              cmz <- max(cmz)*1.5-cmz
              crt <- max(crt)*1.1-crt
              
              # cint <-  log10(int_data[x]) / intlog
              score <- cmz + crt# + cint
              
              if(any(is.na(score))) browser()
              
              return(data.frame(
                from = rep(xi,length(x)),
                to = x + length(mz_ref),
                weight = score,
                mz1=mz_ref[rep(xi,length(x))],
                rt1=rt_ref[rep(xi,length(x))],
                mz2=mz_data[x],
                rt2=rt_data[x]
              ))
            }
            
            # "from","to","weight"
            
          },
          MoreArgs = list(
            mz_data = mz_data,
            rt_data = rt_data,
            mz_ref = mz_ref,
            rt_ref = rt_ref,
            alpha_rt = tol_rt
          ),
          SIMPLIFY = FALSE
        )
      
      
      df_edges <-  do.call(rbind, df_edges)
      ge <-
        make_bipartite_graph(types = c(rep(FALSE, length(mz_ref)), rep(TRUE, length(mz_data))),
                             edges = as.numeric(t(as.matrix(df_edges[,c(1,2)]))),
                             directed = FALSE)
      
      E(ge)$weight <- df_edges[,3]
      
      ###IT WORKS
      mm <- igraph::max_bipartite_match(ge)
      ###We just add the matching part
      res_matching <-  mm$matching[1:length(mz_ref)]
      
      res_matching <-  sapply(res_matching,function(x,shift){
        if(is.na(x)) return(x)
        return(x-shift)
      },shift=length(mz_ref))
      ###Output : a vector giving the postion of mz_ref in the data
      return(res_matching)
    }
  
  matched_features <- matchLCMSsignals(peaktable[,1],
                                     peaktable[,2],
                                     mz,rt,mz_margin,
                                     rt_margin)
  
  matched <- !is.na(matched_features)
  
  ###We select the feature which are matched
  ###we build a feature list.
  lmzr <- split(cbind(mz-mz_margin,mz+mz_margin),f=seq_along(mz))
  lrtr <- split(cbind(rt-rt_margin,rt+rt_margin),f=seq_along(rt))
  
  xraw <- xcmsRaw(path)
  vres <- mapply(lmzr,lrtr,FUN=function(mzr,rtr,xraw){
    reic <- rawEIC(xraw,mzrange=mzr,rtrange=rtr)
    return(list(time=xraw@scantime[reic$scan],intensity=reic$intensity))
  },xraw=xraw)
  

}
