library(xcms)

PATH_COMPOUNDS <- "U:/users/Alexis/data/data_karin/analysis/list_metabo_neg.csv"
compounds <- read.table(PATH_COMPOUNDS,sep=";",header=TRUE,stringsAsFactors = TRUE)

PATH_DATAMATRIX <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res/datamatrices/annotated_peaktable_5c8d1e5ea97664938a069ec6466fbebe_reduced.csv"
dm <- read.table(PATH_DATAMATRIX,sep=";",header=TRUE)


PATH_PEAKTABLES <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res/ADAP/peaktables"
sel_peaktables <- sample(list.files(PATH_PEAKTABLES,full.names=TRUE),15)


path_raw <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/mzML/",samp_name,".mzML",sep="")
xr <- xcmsRaw(path_raw)
pt <- read.table(path_peaktable,sep=",",header=TRUE)


###We check the unique ones.

unalloc_pos <- rep(TRUE,nrow(compounds))


list_group <- list()

for(i in seq_along(unalloc_pos)){
  if(!unalloc_pos[i]) next
  sp <- which(unalloc_pos)
  sel_mass <- which(abs(compounds$M_H_n[sp]-compounds$M_H_n[i])<0.007)
  list_group <- c(list_group,list(sp[sel_mass]))
  unalloc_pos[sp[sel_mass]] <- FALSE
}


###We now extract all the EIC for the associated files 






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
path_peaktable <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/res/ADAP/peaktables/",samp_name,".csv",sep="")
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
    rawEIC(xr,mzrange=c(x[1],x[2]),rtrange=c(rtlims[1],rtlims[2]))
  },xr=xraw)
  return(veics)
}





extractEICsCompounds <- function(xr,peaktable,compounds,tol_mz=0.01,tol_rt=0.3,centered=FALSE){
  source("X:/Documents/dev/script/diverse/match_mz_rt.R")
  
  
  vm <- matchLCMSsignals(mz_data = peaktable[,"mz"],rt_data = peaktable[,"rt"],
                         mz_ref = compounds[,1],rt_ref = compounds[,2],
                         tol_mz = tol_rt,tol_rt = tol_rt)
  
  res <- vector(mode="list",length=nrow(compounds))
  
  found <- which(!is.na(vm))
  
  
  if(length(found)!=0) return(res)
  
  ##The EIC is the peak extended by the peakwidth on the side
  peaktable[]
  
  
  ##We return all the value of the data 
  
  
}







plotCompound




pdf("U:/users/Alexis/data/data_karin/mml_v2/tout.pdf")
for(n_compounds in 1:nrow(compounds)){
  plotCompounds(xr,pt,compounds[n_compounds,1],compounds[n_compounds,"M_H_n"],compounds[n_compounds,"est_RT"])
}
dev.off()