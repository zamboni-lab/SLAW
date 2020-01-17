###Given as input a list of file or a list of ids, extract the associated EICs.
suppressWarnings(suppressMessages(library(DBI)))
suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(BiocParallel)))

ars <- commandArgs(trailingOnly = TRUE)
args <- c(
  "U:/users/Alexis/sandbox/input_features.txt",
  "U:/Alexis/sandbox/input_rawfiles.txt",
  "U:/users/Alexis/examples_lcms_workflow/output/processing_db.sqlite",
  "",
  "",
  "0.02",
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

MARGIN_RT <- 2/60 #(3s)
MARGIN_MZ <- 0.007#MZ tolerance

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


####We check the class of the input
input_signals <- read.table(INPUT_FEATURES,sep=",",header=FALSE,stringsAsFactors = FALSE)


###We read the assocaited data.frame
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
PATH_DATAMATRIX <- dbGetQuery(dbb, "SELECT peaktable FROM PEAKPICKING")[1, 1]
dbDisconnect(dbb)


dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
raw_files <- dbGetQuery(dbb, "SELECT path FROM samples")[, 1]
dbDisconnect(dbb)

dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
peaktables <- dbGetQuery(dbb, "SELECT output_ms FROM processing")[, 1]
dbDisconnect(dbb)


if(FALSE){
  library(stringr)
  PATH_DATAMATRIX <- str_replace(string = PATH_DATAMATRIX,pattern = fixed("/sauer1"),replacement = "U:")
  raw_files <- str_replace(string = raw_files,pattern = fixed("/sauer1"),replacement = "U:")
  peaktables <- str_replace(string = peaktables,pattern = fixed("/sauer1"),replacement = "U:")
}

vmz <- NULL
vrt <- NULL

dm <- read.table(PATH_DATAMATRIX,sep=",",header=TRUE)

if(ncol(input_signals)==1){
  ###Then we consider that the data have been directly removed.
  sel_idx <- input_signals[,1]
  vpp <- log10(apply(dm[,22:ncol(dm)],1,mean,na.rm=TRUE))
  sel_idx <- sample.int(nrow(dm),10,prob = vpp/sum(vpp))
  vmz <- dm[sel_idx,1]
  vrt <- dm[sel_idx,2]
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



titles <- rep("",length(vmz))
if(ncol(input_signals)>=3){
  if("title" %in% colnames(input_signals)){
    titles <- input_signals[,"title"]
  }else{
    titles <- input_signals[,3]
  }
}

supp_titles <- paste("RT: ",sprintf("%0.2f",vrt-MARGIN_RT),"-", sprintf("%0.2f",vrt+MARGIN_RT),
               "\nMZ:",sprintf("%0.4f",vmz-MARGIN_MZ),"-",sprintf("%0.4f",vmz+MARGIN_MZ),sep="")
titles <- paste(titles,supp_titles,sep=" ")




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
  # 
  # ccc <- cbind(peaktable[matched_features[matched],1],
  #       peaktable[matched_features[matched],2],
  #       mz[matched],rt[matched])
  # print(apply(ccc,1,paste,collapse="_"))
  
  ###We select the feature which are matched
  ###we build a feature list.
  rlmzr <- split(cbind(peaktable[matched_features[matched],"mz_min"],
                      peaktable[matched_features[matched],"mz_max"]),f=1:sum(matched))
  rlrtr <- split(cbind((peaktable[matched_features[matched],"rt_min"])*60,
                      (peaktable[matched_features[matched],"rt_max"])*60),f=1:sum(matched))
  
  lmzr <- split(cbind(peaktable[matched_features[matched],"mz_min"]-mz_margin,
                      peaktable[matched_features[matched],"mz_max"]+mz_margin),f=1:sum(matched))
  lrtr <- split(cbind((peaktable[matched_features[matched],"rt_min"]-rt_margin)*60,
                      (peaktable[matched_features[matched],"rt_max"]+rt_margin)*60),f=1:sum(matched))
  xraw <- xcmsRaw(path_xraw)
  reslist <- vector(mode="list",length = length(mz))

  vres <- mapply(lmzr,lrtr,rlmzr,rlrtr,FUN=function(mzr,rtr,rmzr,rrtr,xraw){
    reic <- rawEIC(xraw,mzrange=mzr,rtrange=rtr)
    ###We compute the limit in scan.
    integ <- ifelse(((xraw@scantime[reic$scan]<=rrtr[2]) &
                      (xraw@scantime[reic$scan]>=rrtr[1])),rep(TRUE,length(reic$scan)),
                                                            rep(FALSE,length(reic$scan)))
    
    return(list(time=xraw@scantime[reic$scan],intensity=reic$intensity,peak=integ))
  },MoreArgs = list(xraw=xraw),SIMPLIFY = FALSE)
  
  reslist[matched] <- vres
  
  ### WE return the limits of integration
  
  return(reslist)
}

###We now change the vilsualization of the software by team.



plotPeaks <- function(eics,titles,names_raw){#,name_files,name_compounds

  colors <- rainbow(length(eics))
  
  for(i in seq_along(eics[[1]])){
    rtlim <- sapply(eics,function(x,idx){
      if(is.null(x[[idx]])) return(c(NA_real_,NA_real_,NA_real_))
      c(range(x[[idx]]$time),max(x[[idx]]$intensity))
    },idx=i)
    rtmin <- suppressWarnings(min(rtlim[1,],na.rm=TRUE))
    rtmax <- suppressWarnings(max(rtlim[2,],na.rm=TRUE))
    intmax <- suppressWarnings(max(rtlim[3,],na.rm=TRUE))
    if(is.infinite(rtmin)) next
    
    ####We build the legend vector
    sel_raw <- which(!is.na(rtlim[1,]))
    col_leg <- colors[sel_raw]
    

    plot(0,xlab="Time(min)",ylab="Intensity",xlim=c(rtmin*0.95/60,rtmax*1.05/60),
         ylim=c(0,intmax*1.05),type="n",main=titles[i])
    
    ###We plot all the values
    for(j in seq_along(eics)){
      if(is.null(eics[[j]][[i]])||all(eics[[j]][[i]]$intensity==0)) next
      eic <- eics[[j]][[i]]
      lines(eic$time/60,eic$intensity,lwd=2,col="black")
      lines(eic$time[eic$peak]/60,eic$intensity[eic$peak],lwd=2,col=colors[j])
    }
    legend("topright",legend = names_raw[sel_raw],
           col = colors[sel_raw],lwd = 2,cex = 0.7)
  }
}


bpp <- NULL
if (get_os() == "win") {
  bpp <- SnowParam(workers = NCORES)
} else{
  bpp <- MulticoreParam(workers = min(NCORES, 4))
}

extracted_eics <- bpmapply(raw_files,peaktables,FUN=extract_EIC_raw_file,SIMPLIFY = FALSE,BPPARAM = bpp,
                           MoreArgs = list(mz=vmz,rt=vrt,mz_margin=MARGIN_MZ,rt_margin=MARGIN_RT))

traw_files <- str_split(raw_files,fixed("."),simplify = TRUE)
traw_files <- basename(traw_files[,1])


plotPeaks(extracted_eics,titles = titles,names_raw = traw_files)
# 
# 
# 
# dm[sel_idx[10],c(1,2,22:ncol(dm))]
# 
# xraw <- xcmsRaw(raw_files[2])
# plotEIC(xraw,mzrange=c(dm[sel_idx[3],1]-0.005,dm[sel_idx[3],1]+0.005),
#         rtrange=c(dm[sel_idx[3],2]-0.05,dm[sel_idx[3],2]+0.05)*60)
# rawEIC(xraw,mzrange=c(dm[sel_idx[3],1]-0.005,dm[sel_idx[3],1]+0.005),
#         rtrange=c(dm[sel_idx[3],2]-0.05,dm[sel_idx[3],2]+0.05)*60)$intensity
