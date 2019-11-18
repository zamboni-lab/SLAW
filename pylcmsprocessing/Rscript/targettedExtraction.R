library(xcms)
library(stringr)


###Input the output fuloder of an acquisition

args <- commandArgs(trailingOnly = FALSE)

###Example on the old matrix 

PATH_DATAMATRIX <- args[1]
PATH_DB <- args[2]
PATH_SAMPLES <- args[3]
PATH_COMPOUNDS <- 

TOL_MZ <- as.numeric(args[4])
TOL_RT <- as.numeric(args[5])

###Output a subdataMatrix with the nanotation, aswell asa plot for 1 hundred files.
# * A visulisation of hte targetted peaks on at least 100 samples
# * A  subdata matrix

PATH_COMPOUNDS <- "U:/users/Alexis/data/data_karin/analysis/list_metabo_neg.csv"
compounds <- read.table(PATH_COMPOUNDS,sep=";",header=TRUE,stringsAsFactors = TRUE)

PATH_DATAMATRIX <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res_save/datamatrices/annotated_peaktable_5c8d1e5ea97664938a069ec6466fbebe_reduced.csv"
dm <- read.table(PATH_DATAMATRIX,sep=";",header=TRUE)


samp_name <- "13_SPF_C_c_07"
PATH_PEAKTABLES <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res/ADAP/peaktables"
all_files <- list.files(PATH_PEAKTABLES,full.names=FALSE)
all_samples <- all_files[str_detect(all_files,"[0-9]+_[a-zA-Z]+_[A-E]_[a-m]")]
all_QCs <- all_files[str_detect(all_files,"gX")]



sel_peaktables <- sample(list.files(PATH_PEAKTABLES,full.names=FALSE),100)
path_peaktables <- file.path(PATH_PEAKTABLES,sel_peaktables)
path_samples_peaktables <- file.path(PATH_PEAKTABLES,all_samples)
path_qcs_peaktables <- file.path(PATH_PEAKTABLES,all_QCs)

names_raw <- file.path("U:/users/Alexis/data/data_karin/mml_v2/neg/mzML",str_split(sel_peaktables,"\\.csv",simplify = TRUE)[,1])
names_samples_raw <- file.path("U:/users/Alexis/data/data_karin/mml_v2/neg/mzML",str_split(path_samples_peaktables,"\\.csv",simplify = TRUE)[,1])
names_samples_raw <- paste(names_samples_raw,"mzML",sep=".")
names_qcs_raw <- file.path("U:/users/Alexis/data/data_karin/mml_v2/neg/mzML",str_split(path_qcs_peaktables,"\\.csv",simplify = TRUE)[,1])
names_qcs_raw <- paste(names_qcs_raw,"mzML",sep=".")
# path_raw <- paste(names_raw,".mzML",sep="")
# xr <- xcmsRaw(path_raw[1])
# pt <- read.table(path_peaktable,sep=",",header=TRUE)

RT_TOL <- 0.4
MZ_TOL <- 0.01






###We check the unique ones.
extractUniqueMasses <- function(compounds){
  unalloc_pos <- rep(TRUE,nrow(compounds))
  list_group <- list()
  
  for(i in seq_along(unalloc_pos)){
    if(!unalloc_pos[i]) next
    sp <- which(unalloc_pos)
    sel_mass <- which(abs(compounds$M_H_n[sp]-compounds$M_H_n[i])<0.007)
    list_group <- c(list_group,list(sp[sel_mass]))
    unalloc_pos[sp[sel_mass]] <- FALSE
  }
  return(unalloc_pos)
}



##Extract the heatmap fo the matched ocmpounds
path_peaktables <- path_samples_peaktables[1]
extractPeakTable <- function(path_peaktable,compounds,rt_tol,mz_tol){
  pt <- read.table(path_peaktable,header=TRUE,sep=",")  
  vmm <- matchLCMSsignals(mz_data = pt[,"mz"],rt_data =pt[,"rt"],mz_ref = compounds[,"M_H_n"],rt_ref = compounds[,"est_RT"],
                          tol_rt = rt_tol,tol_mz = mz_tol)
  sel_mol <- !is.na(vmm)
  intval <- rep(NA,nrow(compounds))
  rtval <- rep(NA,nrow(compounds))
  heightval <- rep(NA,nrow(compounds))
  peakwidthval <- rep(NA,nrow(compounds))
  intval[sel_mol] <- pt[vmm[sel_mol],"intensity"]
  heightval[sel_mol] <- pt[vmm[sel_mol],"height"]
  rtval[sel_mol] <- pt[vmm[sel_mol],"rt"]
  peakwidthval[sel_mol] <- pt[vmm[sel_mol],"peakwidth"]
  vsamp <- rep(basename(path_peaktable),nrow(compounds))
  return(data.frame(rt=rtval,intensity=intval,peakwidth=peakwidthval,
                    sample=vsamp,compounds=seq_along(compounds),height=heightval))
}



##For all samples
vam <- sapply(path_samples_peaktables,FUN = extractPeakTable,rt_tol = 0.8,mz_tol=0.01,compounds=compounds,simplify = FALSE)


intMatrix <- do.call(cbind,lapply(vam,FUN = function(x){
  x[,"intensity"]
}))
row.names(intMatrix) <- compounds$name
colnames(intMatrix) <- str_split(basename(names(vam)),pattern = fixed("\\."),simplify = TRUE)[,1]

NA_VAL <- 1000

intMatrix[is.na(intMatrix)] <- NA_VAL
intMatrix <- log10(intMatrix)
library(gplots)
detected <- apply(intMatrix,1,function(x,ref){any(x!=ref)},ref=log10(NA_VAL))

heatmap.2(intMatrix[detected,,drop=FALSE],na.rm = TRUE,col = viridis(begin=0.1,n=100),trace="none",
          dendrogram = "none",Rowv=FALSE)

##Number of detections by sample
num_detect <- apply(intMatrix,1,function(x){sum(x!=log10(NA_VAL))})

boxplot(num_detect)
###Comparison with Karin data
PATH_QUANT_NEG <- "U:/users/Alexis/data/data_karin/mml_v2/neg/quant_karin_neg.csv"
quant_neg <- read.table(PATH_QUANT_NEG,sep=";",header=TRUE,stringsAsFactors = FALSE,na.strings = "")

for(is in 1:ncol(quant_neg)){
  x <- as.character(quant_neg[,is])
  pna <- which(x=="")
  x[pna] <- NA
  x <- as.numeric(x)
  quant_neg[,is] <- x
}


quant_neg <- t(quant_neg)
head(quant_neg[1:5,1:5])
head(intMatrix[1:5,1:5])

vadd <- adist(row.names(quant_neg),row.names(intMatrix))


vm2 <- apply(vadd,1,which.min)
duplicated(vm2)



vam2 <- sapply(path_qcs_peaktables,FUN = extractPeakTable,rt_tol = 0.8,mz_tol=0.01,compounds=compounds,simplify = FALSE)


intMatrix2 <- do.call(cbind,lapply(vam2,FUN = function(x){
  x[,"height"]
}))
row.names(intMatrix2) <- compounds$name
colnames(intMatrix2) <- str_split(basename(names(vam2)),pattern = fixed("\\."),simplify = TRUE)[,1]

NA_VAL <- 1000

intMatrix2[is.na(intMatrix2)] <- NA_VAL
intMatrix2 <- log10(intMatrix2)
library(gplots)
detected <- apply(intMatrix2,1,function(x,ref){any(x!=ref)},ref=log10(NA_VAL))

heatmap.2(intMatrix2[detected,,drop=FALSE],na.rm = TRUE,col = viridis(begin=0.1,n=100),trace="none",
          dendrogram = "none",Rowv=FALSE)

##Number of detections by sample
num_detect <- apply(intMatrix2,1,function(x){sum(x!=log10(NA_VAL))})

boxplot(num_detect)



###Percentage of missing values...
sum(intMatrix!=log10(NA_VAL))*100/length(intMatrix)
###Let s investigate...

vam <- do.call(rbind,)

vmat <- matrix(NA,nrow = nrow(compounds),ncol = length(vam))
vmat[cbind(vam$seq_along.compounds.,vam$sample)]


ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile()


vmm <- matchLCMSsignals(mz_data = dm$mz,rt_data = dm$rt,mz_ref = compounds$M_H_n,rt_ref = compounds$est_RT,tol_mz =  0.01,tol_rt = 1)


###! matching compounds
source("X:/Documents/dev/script/diverse/match_mz_rt.R")

vmm <- matchLCMSsignals(mz_data = dm$mz,rt_data = dm$rt,mz_ref = compounds$M_H_n,rt_ref = compounds$est_RT,tol_mz =  0.01,tol_rt = 1)

vdf <- data.frame(mz=compounds$M_H_n,found=!is.na(vmm))

plotPeaks <- function(peaktable,raw,mz_lim,rt_lim,title=NULL){
  if(class(raw)!="xcmsRaw"){
    raw <- xcmsRaw(raw)
  }
  if(!is.data.frame(peaktable)){
    peaktable <- read.table(peaktable,sep=",",header=TRUE)
  }
  
  
  selp <- which(peaktable$mz>mz_lim[1]&peaktable$mz<mz_lim[2]&
                  peaktable$rt>rt_lim[1]&peaktable$rt<rt_lim[2])
  vcol <- rainbow(length(selp))
  reic <- rawEIC(raw,mzrange=mz_lim,rtrange=rt_lim*60)
  if(all(reic$intensity==0)) return(FALSE)
  teic <- xraw@scantime[reic$scan]/60
  if(is.null(title)){
    title <- paste("EIC for mz: ",paste(sprintf("%0.3f",mz_lim),collapse = "-")," rt: ",
                   paste(sprintf("%0.3f",rt_lim),collapse = "-"),sep="")
  }
  
  plot(teic,reic$intensity,type="l",lwd=2,ylim=c(-1,max(reic$intensity)),col="black",
       xlab="time",ylab="intensity",main=title)
  
  
  if(length(selp)>0){
    for(i in seq_along(selp)){
      ppmin <- which.min(abs(teic-peaktable[selp[i],"rt_min"]))
      ppmax <- which.min(abs(teic-peaktable[selp[i],"rt_max"]))
      tpol <- teic[ppmin:ppmax]
      tpol <- c(tpol,rev(tpol))
      tint <- c(reic$intensity[ppmin:ppmax],rep(0,length(tpol)/2))
      
      ###We select the point which are aligned in the EIC eventually.
      polygon(x=tpol,y = tint,col=vcol[i])
      # segments(x0=pt[selp[i],"rt_min"],x1=pt[selp[i],"rt_max"],y0=0,y1=0,lwd=2,col = vcol[i])
    }
    lines(teic,reic$intensity,lwd=2)
    legend("topleft",legend=sprintf("%0.2f",peaktable[selp,2]),col=vcol,lwd=2)
  }
  return(TRUE)
}

plotCompound <- function(raw,peaktable,name,mz,rt,rt_win=1,mz_win=0.01){
  title <- paste(name,"MZ:",sprintf("%0.3f",mz),"RT:",sprintf("%0.2f",rt))
  plotPeaks(peaktable,raw,mz_lim=c(mz-mz_win,mz+mz_win),rt_lim=c(rt-rt_win,rt+rt_win),title=title)
}


plotCompounds <- function(raw,peaktable,name,mzs,rts,rt_win=1,mz_win=0.01){
  rt_lim <- range(rts)
  rt_lim <- c(rt_lim[1]-rt_win,rt_lim[2]+rt_win)
  mz_lim <- range(mzs)
  mz_lim <- c(mz_lim[1]-mz_win,mz_lim[2]+mz_win)
  title <- paste("compounds with MZ:",sprintf("%0.3f",mzs[1]),"with RT between:",paste(sprintf("%0.2f",rt_lim),collapse="-"),sep="")
  is_plotted <- plotPeaks(peaktable,raw,mz_lim=mz_lim,rt_lim=rt_lim,title=title)
  if(is_plotted){
  colv <- rainbow(length(mzs))

  ##We plot the name and the label for each compounds
  abline(v=rts,col=colv,lwd=2,lty=2)
  legend("topleft",legend = name,lwd = 2,lty=2,col=colv)
  }
}


samp_name <- "12_GF_E_m_13"
path_peaktable <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/res_save/ADAP/peaktables/",samp_name,".csv",sep="")
path_raw <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/mzML/",samp_name,".mzML",sep="")
xr <- xcmsRaw(path_raw)
pt <- read.table(path_peaktable,sep=",",header=TRUE)




###We extract all the EICs for eahc peak tp pllot them simulatneously eventually.
for(i in seq_along(list_group)){
  plotCompounds(xr,pt,compounds$name[list_group[[i]]],
                compounds$M_H_n[list_group[[i]]],compounds$est_RT[list_group[[i]]])
  
}


##extract the WICs given in a matrix
extractEICs <- function(xraw,mzlims,rtlims){
  veics <- apply(cbind(mzlims,rtlims),1,function(x,xr){
    rawEIC(xr,mzrange=c(x[1],x[2]),rtrange=c(x[3],x[4]))
  },xr=xraw)
  return(veics)
}


extractEICsCompounds <- function(xr,peaktable,vm,tol_mz=0.01,tol_rt=0.15){
  res <- vector(mode="list",length=nrow(compounds))
  pfound <- which(!is.na(vm))
  if(length(pfound)==0) return(res)
  mz_lims <- matrix(c(peaktable[vm[pfound],"mz_min"]-0.005,
                    peaktable[vm[pfound],"mz_max"]+0.005),ncol=2)
  rt_lims <- matrix(c((peaktable[vm[pfound],"rt_min"]-peaktable[vm[pfound],"peakwidth"]*0.5)*60,
                                (peaktable[vm[pfound],"rt_max"]+peaktable[vm[pfound],"peakwidth"]*0.5)*60),ncol=2)
  ##The EIC is the peak extended by the peakwidth on the side
  vEICs <- extractEICs(xr,mz_lims,rt_lims)
  for(i in seq_along(vEICs)){
    if(is.null(vEICs[[i]])) next
    vEICs[[i]]$scan <- xr@scantime[vEICs[[i]]$scan]
    vEICs[[i]]$peak <- rep(FALSE,length(vEICs[[i]]$scan))
    pmin <- which(vEICs[[i]]$scan<((peaktable[vm[pfound[i]],"rt_min"]-0.001)*60))
    pmax <- which(vEICs[[i]]$scan>((peaktable[vm[pfound[i]],"rt_max"]+0.001)*60))
    if(length(pmin)==0){
      pmin <- 1
    }else{
      pmin <- pmin[length(pmin)]
    }
    if(length(pmax)==0){
      pmax <- length(vEICs[[i]]$scan)
    }else{
      pmax <- pmax[1]
    }
    
    vEICs[[i]]$peak <- c(pmin,pmax)
    res[[pfound[i]]] <- vEICs[[i]]
  }
  return(res)
}
ref_compounds <- data.frame(mz_neg=compounds$M_H_n,rt=compounds$est_RT,name=compounds$name)
source("X:/Documents/dev/script/diverse/match_mz_rt.R")


peaktables <- sapply(path_peaktables[vidxs],function(x){
  pt <- read.table(x,sep=",",header=TRUE,stringsAsFactors = FALSE)
  return(pt)
},simplify=FALSE)

all_match <- sapply(peaktables,function(x,compounds,tol_rt,tol_mz){

  vm <- matchLCMSsignals(mz_data = pt[,"mz"],rt_data = pt[,"rt"],
                         mz_ref = compounds[,1],rt_ref = compounds[,2],
                         tol_mz = tol_mz,tol_rt = tol_rt)
},simplify=FALSE,compounds=ref_compounds,tol_rt=0.15,tol_mz=MZ_TOL)


vidxs <- sample(seq_along(path_peaktables),15)


exEICs <- function(x,y,compounds,mz_tol,rt_tol){
  pt <- read.table(x,sep=",",header=TRUE)
  xr <- xcmsRaw(y)
  vm <- matchLCMSsignals(mz_data = pt[,"mz"],rt_data = pt[,"rt"],
                         mz_ref = compounds[,1],rt_ref = compounds[,2],
                         tol_mz = mz_tol,tol_rt =rt_tol)
  vEICs <- extractEICsCompounds(xr,pt,vm,ref_compounds)
}

cEICs <- mapply(path_peaktables[vidxs[1:5]],path_raw[vidxs[1:5]],FUN = exEICs,
                MoreArgs = list(compounds=ref_compounds,mz_tol=MZ_TOL,rt_tol=0.1),SIMPLIFY = FALSE)
# cEICt <- exEICs(path_peaktables[vidxs[1]],path_raw[vidxs[1]],ref_compounds,MZ_TOL,0.1)




plotStackedEICs <- function(cEICs,compounds,sampnames){
  ##We extract all the idx with good vlaue
  
  ###We gives an lwd
  #10 Colors and 5 strokes evnetually and that s it.
  
  vcol <- scales:::hue_pal()(min(length(sampnames),10))
  
  
  for(i in 1:nrow(compounds)){
    cComp <- sapply(cEICs,function(x,idx){
      if(is.null(x[[idx]])){
        return(NULL)
      }
      return(x[[idx]])
    },idx=i,simplify=FALSE)
    
    sel <- which(sapply(cComp,function(x){!is.null(x)}))
    if(length(sel)==0) next
    int_min <- 0
    int_max <- max(sapply(cComp[sel],function(x){max(x[[2]])}))
    rt_min <- min(sapply(cComp[sel],function(x){min(x[[1]])}))
    rt_max <- max(sapply(cComp[sel],function(x){max(x[[1]])}))
    margin <- 0.1*(rt_max-rt_min)
    # if(compounds[i,"name"]=="Glucose 6-phosphate") browser()
    plot(-1000,-1000,xlab="RT",ylab="intensity",xlim=c(rt_min-margin,rt_max+margin)/60,ylim=c(int_min,int_max),
         main=paste(compounds[i,"name"]," (MZ:",sprintf("%0.3f",compounds[i,1]),
                    " RT:",sprintf("%0.2f",compounds[i,"rt"]),")",sep=""))
    for(ss in sel){
      # pcol <- rep("#000000",length(cComp[[ss]][[3]]))
      # pcol[cComp[[ss]][[3]]] <- vcol[ss]
      lines(cComp[[ss]][[1]]/60 ,cComp[[ss]][[2]],type="l",col=vcol[ss],lwd=2)
      lines(cComp[[ss]][[1]][1:cComp[[ss]][[3]][1]]/60,cComp[[ss]][[2]][1:cComp[[ss]][[3]][1]],
            col="black",lwd=2)
      lines(cComp[[ss]][[1]][cComp[[ss]][[3]][2]:length(cComp[[ss]][[1]])]/60,
            cComp[[ss]][[2]][cComp[[ss]][[3]][2]:length(cComp[[ss]][[1]])],
            col="black",lwd=2)
    }
    legend("topleft",legend=sampnames[sel],lwd=2,col=vcol[sel],cex=0.7)
  }
}

plotStackedEICs(cEICs,ref_compounds,basename(path_peaktables[vidxs]))



###Extract the usbmatrix with th



###We extract the table of table for the gXMatrix 

##pos_pe
list_peaktables <- list.files("U:/users/Alexis/data/data_karin/mml_v2/neg/res/ADAP/peaktables",full.names=TRUE)

###We parse lal the peaktables and extract hte measured RT amsses en intensity eventually.

m





pdf("U:/users/Alexis/data/data_karin/mml_v2/tout.pdf")
for(n_compounds in 1:nrow(compounds)){
  plotCompounds(xr,pt,compounds[n_compounds,1],compounds[n_compounds,"M_H_n"],compounds[n_compounds,"est_RT"])
}
dev.off()