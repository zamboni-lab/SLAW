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

PATH_DATAMATRIX <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res_thin/datamatrices/datamatrix_b01904e0a545f19857549535352af3b9.csv"
dm <- read.table(PATH_DATAMATRIX,sep=",",header=TRUE)

"U:/users/Alexis/transfer/for_karin/raw_datamatrix_no_blank.csv"



samp_name <- "13_SPF_C_c_07"
PATH_PEAKTABLES <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res/ADAP/peaktables"
all_files <- list.files(PATH_PEAKTABLES,full.names=FALSE)
all_samples <- all_files[str_detect(all_files,"[0-9]+_[a-zA-Z]+_[A-E]_[a-m]")]
all_QCs <- all_files[str_detect(all_files,"gX")]


sel_peaktables <- sample(list.files(PATH_PEAKTABLES,full.names=FALSE),100)
path_peaktables <- file.path(PATH_PEAKTABLES,sel_peaktables)
path_samples_peaktables <- file.path(PATH_PEAKTABLES,all_samples)

ninj_name <- paste(str_split(path_samples_peaktables,pattern=fixed("."),simplify = TRUE)[,1],"_newinj.csv",sep="")
path_samples_peaktables[file.exists(ninj_name)] <- ninj_name[file.exists(ninj_name)]
path_samples_peaktables <- path_samples_peaktables[!duplicated(path_samples_peaktables)]

##If a file with new injection exists we replace it.



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
  return(list_group)
}

##Extract the heatmap fo the matched ocmpounds
path_peaktables <- path_samples_peaktables[1]
extractPeakTable <- function(path_peaktable,compounds,rt_tol,mz_tol){
  pt <- read.table(path_peaktable,header=TRUE,sep=",")  
  source("X:/Documents/dev/script/diverse/match_mz_rt.R")
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

###We ocnstruct the quantiofication table for xcms
# stst_b <- which(str_detect(line_xcmss,"c\\("))
# stst_e <- which(str_detect(line_xcmss,"\\)"))
# sel_pos_1 <- which(stst_b!=stst_e)
# 
# to_change <- mapply(stst_b[sel_pos_1],stst_e[sel_pos_1],FUN=function(x,y,vline){
#   str_replace(paste(line_xcmss[x:y],collapse=""),"c\\(.+\\)","NA")
# },MoreArgs = list(vline=line_xcmss))
# 
# 
# idx_to_rm <- do.call(c,mapply(stst_b[sel_pos_1],stst_e[sel_pos_1],FUN=function(x,y){
#   (x+1):y
# }))
# 
# line_xcmss[stst_b[sel_pos_1]] <- to_change
# line_xcmss <- line_xcmss[-idx_to_rm]
# 
# writeLines(line_xcmss,con = file("U:/users/Alexis/data/data_karin/mml_v2/neg/res_mml2_xcms_corrected.csv","w"))
line_xcmss <- readLines(con = file(PATH_DM_XCMS,"r"))
line_xcmss <- str_replace(line_xcmss,"c\\(.+\\)","NA")
writeLines(libe_xcmss,con = file("U:/users/Alexis/data/data_karin/mml_v2/neg/res_mml2_xcms_corrected.csv","w"))

PATH_DM_XCMS <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res_mml2_xcms_corrected.csv"
PATH_DM_XCMS <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res_mml2_xcms.csv"
dm_xcms <- read.table(PATH_DM_XCMS,sep=',',header=TRUE)


vmxmcs <- matchLCMSsignals(mz_data = dm_xcms[,"mzmed"],rt_data =dm_xcms[,"rtmed"]/60,mz_ref = compounds[,"M_H_n"],rt_ref = compounds[,"est_RT"],
                        tol_rt = 0.1,tol_mz = 0.01)

vfm <- xcms:::fastMatch(compounds[,"M_H_n"],dm_xcms[,"mzmed"],tol = 0.005)


dm_xcms[vfm[[122]],c("mzmed","rtmed","npeaks")]
compounds[122,]



xcms_matrix_r <- mapply(vfm,as.list(compounds$est_RT),FUN=function(x,rt,dm,tol_rt){
  if(is.null(x)) return(rep(NA,470-10+1))
  sel_pos <- which(abs(dm[x,"rtmed"]/60-rt)<tol_rt)
  if(length(sel_pos)==0) return(rep(NA,470-10+1))
  apply(dm[x[sel_pos],10:470,drop=FALSE],2,max,na.rm=TRUE)
},MoreArgs = list(dm=dm_xcms,tol_rt=0.15),SIMPLIFY = FALSE)

xcms_matrix_r <- do.call(rbind,xcms_matrix_r)
xcms_matrix_r[is.infinite(xcms_matrix_r)] <- NA


row.names(xcms_matrix_r) <- compounds$name

###WE tryt to match the signal in the xcms matrix
# sub_matrix_xcms <- sub_matrix_xcms[,10:470]
cnames <- sapply(colnames(xcms_matrix_r),function(x){
  str_sub(x,2,-6)
})

colnames(xcms_matrix_r) <- cnames
removed_pos <- which(endsWith(cnames,"newinj"))+1

xcms_matrix_r <- xcms_matrix_r[,-removed_pos]
cnames <- colnames(xcms_matrix_r)
wrong_pos <- which(str_detect(cnames,"newinj"))
cnames[wrong_pos] <- str_split(cnames[wrong_pos],fixed("_newinj"),simplify = TRUE)[,1]
colnames(xcms_matrix_r) <- cnames
svmm <- match(colnames(quant_neg),colnames(xcms_matrix_r))
xcms_matrix_r <- xcms_matrix_r[,svmm]

vadd <- adist(row.names(quant_neg),row.names(xcms_matrix_r))
vm2 <- apply(vadd,1,which.min)
xcms_matrix_r <- xcms_matrix_r[vm2,]



dim(xcms_matrix_r)

sum(!is.na(vmxmcs))
length(vmxmcs)
head(vmxmcs)
dim(intMatrix)



###We build a ocmbined matrix 
# sub_matrix_xcms <- dm_xcms[vmxmcs[sel_comp_xcms],]


# row.names(sub_matrix_xcms) <- compounds$name[sel_comp_xcms]

###WE tryt to match the signal in the xcms matrix
# sub_matrix_xcms <- sub_matrix_xcms[,10:470]
# cnames <- sapply(colnames(sub_matrix_xcms),function(x){
#   str_sub(x,2,-6)
# })
# 
# colnames(sub_matrix_xcms) <- cnames
# removed_pos <- which(endsWith(cnames,"newinj"))+1
# 
# sub_matrix_xcms <- sub_matrix_xcms[,-removed_pos]
# cnames <- colnames(sub_matrix_xcms)
# wrong_pos <- which(str_detect(cnames,"newinj"))
# cnames[wrong_pos] <- str_split(cnames[wrong_pos],fixed("_newinj"),simplify = TRUE)[,1]
# colnames(sub_matrix_xcms) <- cnames
# svmm <- match(colnames(quant_neg),colnames(sub_matrix_xcms))
# sub_matrix_xcms <- sub_matrix_xcms[,svmm]
# 
# vadd <- adist(row.names(quant_neg),row.names(sub_matrix_xcms))
# vm2 <- apply(vadd,1,which.min)
# sub_matrix_xcms <- sub_matrix_xcms[vm2,]

sub_matrix_xcms <- xcms_matrix_r

write.table(sub_matrix_xcms,file="U:/users/Alexis/data/data_karin/mml_v2/neg/in_house_quant_xcms.csv",row.names = TRUE,sep=",")






####Building a detection table at the peaktable level.
vam <- bplapply(as.list(path_samples_peaktables),FUN = extractPeakTable,rt_tol = 0.2,mz_tol=0.01,compounds=compounds)

intMatrix <- do.call(cbind,lapply(vam,FUN = function(x){
  x[,"intensity"]
}))

peakwidthMatrix <- do.call(cbind,lapply(vam,FUN = function(x){
  x[,"peakwidth"]
}))

row.names(intMatrix) <- compounds$name
colnames(intMatrix) <- str_split(basename(path_samples_peaktables),pattern = fixed("."),simplify = TRUE)[,1]

cnames <- colnames(intMatrix)
wrong_pos <- which(str_detect(cnames,"newinj"))
cnames[wrong_pos] <- str_split(cnames[wrong_pos],fixed("_newinj"),simplify = TRUE)[,1]
colnames(intMatrix) <- cnames
svmm <- match(colnames(quant_neg),colnames(intMatrix))



PATH_QUANT_NEG <- "U:/users/Alexis/data/data_karin/mml_v2/neg/quant_karin_neg.csv"
quant_neg <- read.table(PATH_QUANT_NEG,sep=";",header=TRUE,stringsAsFactors = FALSE,na.strings = "",row.names = 1)


for(is in 1:ncol(quant_neg)){
  x <- as.character(quant_neg[,is])
  pna <- which(x=="")
  x[pna] <- NA
  x <- as.numeric(x)
  quant_neg[,is] <- x
}
quant_neg <- t(quant_neg)


quant_neg[quant_neg<min(intMatrix,na.rm = TRUE)] <- NA
cnames <- colnames(intMatrix)
wrong_pos <- which(str_detect(cnames,"newinj"))
cnames[wrong_pos] <- str_split(cnames[wrong_pos],fixed("_newinj"),simplify = TRUE)[,1]
colnames(intMatrix) <- cnames
svmm <- match(colnames(quant_neg),colnames(intMatrix))

# svmm <- match(paste(colnames(quant_neg),"_newinj",sep=""),colnames(intMatrix))
intMatrix <- intMatrix[,svmm]

write.table(intMatrix,file="U:/users/Alexis/data/data_karin/mml_v2/neg/in_house_quant_thin_2.csv",row.names = TRUE,sep=",")

###We noew reorder the rows

vadd <- adist(row.names(quant_neg),row.names(intMatrix))
vm2 <- apply(vadd,1,which.min)
intMatrix <- intMatrix[vm2,]

###We remove all the rest eventually.

lintMatrix <- intMatrix
NA_VAL <- 10
lintMatrix[is.na(lintMatrix)] <- NA_VAL
lintMatrix <- log10(lintMatrix)

lquant <- quant_neg
NA_VAL <- 10
lquant[is.na(lquant)] <- NA_VAL
lquant <- log10(lquant)


lxcms <- as.matrix(xcms_matrix_r)
NA_VAL <- 10
lxcms[is.na(lxcms)] <- NA_VAL
lxcms <- log10(lxcms)

library(gplots)
library(viridis)

order_int <- order(apply(lquant,1,sum,na.rm=TRUE),decreasing = TRUE)
o_sample <- order(apply(lquant,2,sum,na.rm=TRUE),decreasing = TRUE)
heatmap.2(lquant[order_int,o_sample],na.rm = TRUE,col = viridis(begin=0.1,n=100),trace="none",
          dendrogram = "none",Rowv=FALSE,Colv = FALSE)
heatmap.2(lintMatrix[order_int,o_sample],na.rm = TRUE,col = viridis(begin=0.1,n=100),trace="none",
          dendrogram = "none",Rowv=FALSE,Colv = FALSE)
heatmap.2(lxcms[order_int,o_sample],na.rm = TRUE,col = viridis(begin=0.1,n=100),trace="none",
          dendrogram = "none",Rowv=FALSE,Colv = FALSE)
library(gplots)
detected <- apply(intMatrix,1,function(x,ref){any(x!=ref)},ref=log10(NA_VAL))

heatmap.2(intMatrix[detected,,drop=FALSE],na.rm = TRUE,col = viridis(begin=0.1,n=100),trace="none",
          dendrogram = "none",Rowv=FALSE)

###We check the nuber of detectection in both dataset
bintMatrix <- (!is.na(intMatrix))*2
bintQuant <- (!is.na(quant_neg))
binPres <- bintMatrix+bintQuant

###log10 intensity

sintM <- which(!is.na(intMatrix),arr.ind = TRUE)
sintQ <- which(!is.na(quant_neg),arr.ind = TRUE)


df_infos<- data.frame(sample=c(sintM[,2],sintQ[,2]),
                      compounds=c(sintM[,1],sintQ[,1]),
                      int=c(lintMatrix[sintM],lquant[sintQ]),
                      type=as.factor(c(rep("LC-MS workflow",nrow(sintM)),
                                       rep("Quant-Karin",nrow(sintQ)))))
onlyM <- which(binPres==2,arr.ind = TRUE)
onlyQ <- which(binPres==1,arr.ind = TRUE)
both <- which(binPres==3,arr.ind = TRUE)

df_infos_both <- data.frame(sample=c(onlyM[,2],onlyQ[,2],both[,2]),
                            compounds=c(onlyM[,1],onlyQ[,1],both[,1]),
                            int=c(lintMatrix[onlyM],lquant[onlyQ],lquant[both]),
                            detection=as.factor(c(rep("cWorfklow only",nrow(onlyM)),
                                             rep("Targetted only",nrow(onlyQ)),
                                             rep("Both",nrow(both)))))


library(ggplot2)
ggp <- ggplot(df_infos,aes(x=int,fill=type),alpha=0.1)+geom_histogram(position="stack",bins=50)+xlab("Log10(count)")+ylab("Number of features")+
  theme(legend.title = element_text(size=12),legend.text=element_text(size=12))+labs(fill="Type of quantification")


ggp <- ggplot(df_infos,aes(x=int,y=..scaled..,fill=type),alpha=0.1)+geom_density(alpha=0.5)+xlab("Log10(Area)")+ylab("Number of features")+
  theme(legend.title = element_text(size=12),legend.text=element_text(size=12))+labs(fill="Type of quantification")
plot(ggp)

png("U:/users/Alexis/presentation/Group_meeting_25_11/density_karin_quant_thin.png")
plot(ggp)
dev.off()

ggp <- ggplot(df_infos_both,aes(x=int,y=..count..,fill=detection),
              alpha=0.1)+geom_density(alpha=0.5,bw=0.1)+xlab("Log10(Area)")+ylab("Number of features")+
  theme(legend.title = element_text(size=12),legend.text=element_text(size=12))+labs(fill="Detected by")
png("U:/users/Alexis/presentation/Group_meeting_25_11/density_detection_quant_thin.png")
plot(ggp)
dev.off()


###XCMS_INFOS <- 

bintMatrixX <- (!is.na(xcms_matrix_r))*2
binPresX <- bintMatrixX+bintQuant
###log10 intensity
sintMX <- which(!is.na(xcms_matrix_r),arr.ind = TRUE)



df_infos<- data.frame(sample=c(sintMX[,2],sintQ[,2]),
                      compounds=c(sintMX[,1],sintQ[,1]),
                      int=c(lxcms[sintMX],lquant[sintQ]),
                      type=as.factor(c(rep("XCMS",nrow(sintMX)),
                                       rep("Quant-Karin",nrow(sintQ)))))
onlyMX <- which(binPresX==2,arr.ind = TRUE)
# onlyQ <- which(binPres==1,arr.ind = TRUE)
bothX <- which(binPresX==3,arr.ind = TRUE)

###For the summary density
bothLX <- which(binPresX==3&binPres==3,arr.ind = TRUE)
xcmsOnly <- which(binPresX==3&binPres==1,arr.ind = TRUE)
wfOnly <- which(binPres==3&binPresX==1,arr.ind = TRUE)
##Missed by both 
missedBoth <- which(binPresX==1&binPres==1&bintQuant,arr.ind = TRUE)

df_infos_summary <- data.frame(sample=c(bothLX[,2],xcmsOnly[,2],wfOnly[,2],missedBoth[,2]),
                               compounds=c(bothLX[,1],xcmsOnly[,1],wfOnly[,1],missedBoth[,1]),
                               int=c(lxcms[bothLX],lxcms[xcmsOnly],lintMatrix[wfOnly],lquant[missedBoth]),
                               detection=as.factor(c(rep("All methods",nrow(bothLX)),
                                                     rep("XCMS and\n targetted",nrow(xcmsOnly)),
                                                     rep("In-house and\n targetted",nrow(wfOnly)),
                                                     rep("targetted only",nrow(missedBoth)))))


ggp <- ggplot(df_infos_summary,aes(x=int,y=..count..,fill=detection),
              alpha=0.1)+geom_density(alpha=0.5,bw=0.1)+xlab("Log10(Area)")+ylab("Number of features")+
  theme(legend.title = element_text(size=15),legend.text=element_text(size=15))+labs(fill="Detected by")
png("U:/users/Alexis/presentation/Group_meeting_25_11/density_detection_summary.png")
plot(ggp)
dev.off()

df_infos_bothX <- data.frame(sample=c(onlyMX[,2],onlyQ[,2],bothX[,2]),
                            compounds=c(onlyMX[,1],onlyQ[,1],bothX[,1]),
                            int=c(lxcms[onlyMX],lquant[onlyQ],lquant[bothX]),
                            detection=as.factor(c(rep("cXCMS",nrow(onlyMX)),
                                                  rep("Targetted only",nrow(onlyQ)),
                                                  rep("Both",nrow(bothX)))))




library(ggplot2)
ggp <- ggplot(df_infos_bothX,aes(x=int,fill=detection),alpha=0.1)+geom_histogram(position="stack",bins=50)+xlab("Log10(count)")+ylab("Number of features")+
  theme(legend.title = element_text(size=12),legend.text=element_text(size=12))+labs(fill="Type of quantification")


ggp <- ggplot(df_infos_bothX,aes(x=int,y=..scaled..,fill=detection),alpha=0.1)+geom_density(alpha=0.5)+xlab("Log10(Area)")+ylab("Number of features")+
  theme(legend.title = element_text(size=12),legend.text=element_text(size=12))+labs(fill="Type of quantification")
plot(ggp)

png("U:/users/Alexis/presentation/Group_meeting_25_11/density_karin_quant_thin_xcms.png")
plot(ggp)
dev.off()

ggp <- ggplot(df_infos_bothX,aes(x=int,y=..count..,fill=detection),
              alpha=0.1)+geom_density(alpha=0.5,bw=0.1)+xlab("Log10(Area)")+ylab("Number of features")+
  theme(legend.title = element_text(size=12),legend.text=element_text(size=12))+labs(fill="Detected by")
png("U:/users/Alexis/presentation/Group_meeting_25_11/density_detection_quant_thin_xcms.png")
plot(ggp)
dev.off()



##Errro standard data 
###CHekcing the undetected onw
error_mat <- matrix(NA_real_,nrow=nrow(intMatrix),ncol=ncol(intMatrix))
error_mat[both] <- abs(intMatrix[both]-quant_neg[both])/quant_neg[both]

sel_error <- which(!is.na(error_mat),arr.ind=TRUE)

vbin <- .bincode(lquant[sel_error],breaks = 1:10)
tdff <- data.frame(rel_error=error_mat[sel_error],log_error=log10(error_mat[sel_error]),
                   compounds=row.names(intMatrix)[sel_error[,1]],int = lquant[sel_error],bin_int=vbin)
tdff$log_error[tdff$log_error<= -2] <- -2
ggp <- ggplot(tdff,aes(y=log_error,x=as.factor(bin_int)))+geom_boxplot()
plot(ggp)

###Eroro xcms 
error_matX <- matrix(NA_real_,nrow=nrow(xcms_matrix_r),ncol=ncol(xcms_matrix_r))
error_matX[bothX] <- abs(xcms_matrix_r[bothX]-quant_neg[bothX])/quant_neg[bothX]

sel_errorX <- which(!is.na(error_matX),arr.ind=TRUE)

vbinX <- .bincode(lquant[sel_errorX],breaks = 1:10)
tdffX <- data.frame(rel_error=error_matX[sel_errorX],log_error=log10(error_matX[sel_errorX]),
                   compounds=row.names(xcms_matrix_r)[sel_errorX[,1]],int = lquant[sel_errorX],bin_int=vbinX)
tdffX$log_error[tdffX$log_error<= -2] <- -2

ggp <- ggplot(tdffX,aes(y=log_error,x=as.factor(bin_int)))+geom_boxplot()
plot(ggp)



fun_count <- function(x){
  med <- median(10**x)
  if(med==1){
    med <- 0
  }
  return(data.frame(label=paste("count:",length(x),"\nmedian:",sprintf("%0.1f",med),
                                "\nmean:",sprintf("%0.1f",mean(10**x)),sep=""),y=max(x)+0.4))
}

fun_median <- function(x){
  return(data.frame(label=sprintf("%0.2f",median(x)),y=median(x)+0.2))
}

ggp <- ggplot(tdff,aes(y=rel_error*100,x=as.factor(bin_int),col=as.factor(bin_int)))+geom_boxplot()+
  xlab("log10(Area quant)")+ylab("Relative error of integration (%)")+scale_color_viridis_d(end = 0.85)+
  stat_summary(fun.data = fun_count, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75),col="black")+guides(color = FALSE)

png("U:/users/Alexis/presentation/Group_meeting_25_11/error_quant.png")
plot(ggp)
dev.off()

ggp <- ggplot(tdffX,aes(y=rel_error*100,x=as.factor(bin_int),col=as.factor(bin_int)))+geom_boxplot()+
  xlab("log10(Area quant)")+ylab("Relative error of integration (%)")+scale_color_viridis_d(end = 0.85)+
  stat_summary(fun.data = fun_count, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75),col="black")+guides(color = FALSE)

png("U:/users/Alexis/presentation/Group_meeting_25_11/error_quant.png")
plot(ggp)
dev.off()


ggp <- ggplot(tdff,aes(y=log_error+2,x=as.factor(bin_int),col=as.factor(bin_int)))+geom_boxplot()+
  xlab("log10(Area quant)")+ylab("log(10) Relative error of integration (%)")+scale_color_viridis_d(end = 0.85)+
  stat_summary(fun.data = fun_count, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75),col="black")+guides(color = FALSE)+ylim(c(0,4.5))


# ggp <- ggplot(tdff,aes(y=log_error+2,x=as.factor(bin_int),col=as.factor(bin_int)))+geom_violin()+geom_boxplot(width=0.1)+
#   xlab("log10(Area quant)")+ylab("log(10) Relative error of integration (%)")+scale_color_viridis_d(end = 0.85)+
#   stat_summary(fun.data = fun_count, geom = "text", fun.y = median, 
#                position = position_dodge(width = 0.75),col="black")+guides(color = FALSE)

png("U:/users/Alexis/presentation/Group_meeting_25_11/log_error_quant_workflow_boxplot.png")
plot(ggp)
dev.off()

sel_errorXB <- which(!is.na(error_matX)&!is.na(error_mat),arr.ind=TRUE)
vbinXB <- .bincode(lquant[sel_errorXB],breaks = 1:10)
tdffXB <- data.frame(rel_error=error_matX[sel_errorXB],log_error=log10(error_matX[sel_errorXB]),
                    compounds=row.names(xcms_matrix_r)[sel_errorXB[,1]],int = lquant[sel_errorXB],bin_int=vbinXB)
tdffB <- data.frame(rel_error=error_mat[sel_errorXB],log_error=log10(error_mat[sel_errorXB]),
                     compounds=row.names(intMatrix)[sel_errorXB[,1]],int = lquant[sel_errorXB],bin_int=vbinXB)
tdffXB$log_error[tdffXB$log_error<= -2] <- -2
tdffB$log_error[tdffB$log_error<= -2] <- -2
library(scales)
ggp <- ggplot(tdffXB,aes(y=log_error+2,x=as.factor(bin_int),col=as.factor(bin_int)))+geom_boxplot()+
  xlab("log10(Area quant)")+ylab("log(10) Relative error of integration (%)")+scale_color_viridis_d(end = 0.85)+
  stat_summary(fun.data = fun_count, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75),col="black")+guides(color = FALSE)+ylim(c(0,4.5))

png("U:/users/Alexis/presentation/Group_meeting_25_11/log_error_quant_xcms_boxplot_both.png")
plot(ggp)
dev.off()

ggp <- ggplot(tdffB,aes(y=log_error+2,x=as.factor(bin_int),col=as.factor(bin_int)))+geom_boxplot()+
  xlab("log10(Area quant)")+ylab("log(10) Relative error of integration (%)")+scale_color_viridis_d(end = 0.85)+
  stat_summary(fun.data = fun_count, geom = "text", fun.y = median,position = position_dodge(width = 0.75),col="black")+guides(color = FALSE)+ylim(c(0,4.5))
png("U:/users/Alexis/presentation/Group_meeting_25_11/log_error_quant_workflow_boxplot_both.png")
plot(ggp)
dev.off()







sel_comp <- which.max(apply(error_mat,1,mean,na.rm=TRUE))

plot(lintMatrix[sel_comp,],col="blue",pch=2)
points(lquant[sel_comp,],col="red")


sel_sample <- which.max(error_mat[sel_comp,])
rownames(intMatrix)[sel_comp]
colnames(intMatrix)[sel_sample]

intMatrix[sel_comp,sel_sample]
intMatrix[sel_comp,sel_sample]



###We plot the complex region
z1 <- c('Ursodeoxycholic acid','Hyodeoxycholic acid',
        'Chenodeoxycholic acid','Deoxycholic acid')
ssamp <- 1:ncol(intMatrix)




sel_col_z1 <- apply(lquant[match(z1,rownames(intMatrix)),ssamp],2,FUN=function(x){any(x!=1)})

os1 <- order(apply(10**lquant[match(z1,rownames(intMatrix)),sel_col_z1],2,sum,na.rm=TRUE),decreasing=TRUE)

png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_4_w.png")
heatmap.2(t(lintMatrix[match(z1,rownames(intMatrix)),sel_col_z1][,os1]),trace = "none",col = viridis,Colv=FALSE,Rowv = FALSE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE,key=FALSE)
dev.off()
png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_4_quant.png")
heatmap.2(t(lquant[match(z1,rownames(intMatrix)),sel_col_z1][,os1]),trace = "none",col = viridis,Colv=FALSE,Rowv = FALSE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE,key=FALSE)
dev.off()


png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_4_x.png")
heatmap.2(t(lxcms[match(z1,rownames(xcms_matrix_r)),sel_col_z1][,os1]),trace = "none",col = viridis,Colv=FALSE,Rowv = FALSE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE,key=FALSE)
dev.off()


z2 <- c('a-Muricholic acid','Hyocholic acid','Cholic acid')
sel_col_z2 <- apply(lquant[match(z2,rownames(intMatrix)),ssamp],2,FUN=function(x){any(x!=1)})
os2 <- order(apply(10**lquant[match(z2,rownames(intMatrix)),sel_col_z2],2,sum,na.rm=TRUE),decreasing=TRUE)

png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_3_w.png")
heatmap.2(t(lintMatrix[match(z2,rownames(intMatrix)),sel_col_z2][,os2]),trace = "none",col = viridis,Colv=FALSE,Rowv = FALSE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE,key=FALSE)
dev.off()
png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_3_quant.png")
heatmap.2(t(lquant[match(z2,rownames(intMatrix)),sel_col_z2][,os2]),trace = "none",col = viridis,Colv=FALSE,Rowv = FALSE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE,key=FALSE)
dev.off()


png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_3_x.png")
heatmap.2(t(lxcms[match(z2,rownames(xcms_matrix_r)),sel_col_z2][,os2]),trace = "none",col = viridis,Colv=FALSE,Rowv = FALSE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE,key=FALSE)
dev.off()


match(z2,rownames(intMatrix))
ssamp <- sample.int(nrow(intMatrix),100)
png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_3_w.png")
heatmap.2(t(lintMatrix[match(z2,rownames(intMatrix)),ssamp]),trace = "none",col = viridis,Colv=FALSE,Rowv = TRUE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE)
dev.off()

png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_3_quant_w.png")
heatmap.2(t(lquant[match(z2,rownames(intMatrix)),ssamp]),trace = "none",col = viridis,Colv=FALSE,Rowv = TRUE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE)
dev.off()

sel_col_z_1_x <- apply(lquant[match(z1,rownames(xcms_matrix_r)),ssamp],2,FUN=function(x){any(x!=1)})


png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_4_quant_x.png")
heatmap.2(t(lquant[match(z1,rownames(lxcms)),sel_col_z_1_x]),trace = "none",col = viridis,Colv=FALSE,Rowv = FALSE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE)
dev.off()

z2 <- c('a-Muricholic acid','Hyocholic acid','Cholic acid')
match(z2,rownames(xcms_matrix_r))
ssamp <- sample.int(nrow(xcms_matrix_r),100)
png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_3_x.png")
heatmap.2(t(lxcms[match(z2,rownames(xcms_matrix_r)),ssamp]),trace = "none",col = viridis,Colv=FALSE,Rowv = TRUE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE)
dev.off()

png("U:/users/Alexis/presentation/Group_meeting_25_11/iso_3_quant_x.png")
heatmap.2(t(lquant[match(z2,rownames(xcms_matrix_r)),ssamp]),trace = "none",col = viridis,Colv=FALSE,Rowv = TRUE,
          cexRow=0.2,dendrogram="none",lmat= rbind(c(0,3),c(0,1),c(4,1),c(0,1),c(0,2)),lhei=c(0.6,1.5,1.5,1.5,0.6),
          labCol=FALSE,labRow=FALSE)
dev.off()




heatmap.2(lquant[match(z2,rownames(intMatrix)),ssamp],trace = "none",col = viridis,Colv = TRUE,Rowv = FALSE)
plot(peakwidthMatrix[match(z2[2],rownames(intMatrix)),ssamp],lintMatrix[match(z2[2],rownames(intMatrix)),ssamp])
###! matching compounds
source("X:/Documents/dev/script/diverse/match_mz_rt.R")

vmm <- matchLCMSsignals(mz_data = dm$mz,rt_data = dm$rt,mz_ref = compounds$M_H_n,rt_ref = compounds$est_RT,tol_mz =  0.01,tol_rt = 1)

vdf <- data.frame(mz=compounds$M_H_n,found=!is.na(vmm))

plotPeaks <- function(peaktable,raw,mz_lim,rt_lim,title=NULL,addLeg=TRUE){
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
    area_sum <- numeric(length(selp))
    area_normalized <- numeric(length(selp))
    for(i in seq_along(selp)){
      ppmin <- which.min(abs(teic-peaktable[selp[i],"rt_min"]))-1
      ppmax <- which.min(abs(teic-peaktable[selp[i],"rt_max"]))
      tpol <- teic[ppmin:ppmax]
      tpol <- c(tpol,rev(tpol))
      tint <- c(reic$intensity[ppmin:ppmax],rep(0,length(tpol)/2))
      area_sum[i] <- sum(reic$intensity[(ppmin+1):ppmax])
      area_normalized[i] <- trapz(teic[(ppmin+1):ppmax]*60,y = reic$intensity[(ppmin+1):ppmax])
      ###We select the point which are aligned in the EIC eventually.
      polygon(x=tpol,y = tint,col=vcol[i])
      
      ###We calculate the two kind of integration, one while summing everythig
      # segments(x0=pt[selp[i],"rt_min"],x1=pt[selp[i],"rt_max"],y0=0,y1=0,lwd=2,col = vcol[i])
    }
    lines(teic,reic$intensity,lwd=2)
    if(addLeg)
    legend("topleft",legend=paste("RT:",sprintf("%0.2f",peaktable[selp,2])," SUM:",
                                  sprintf("%0.2f",area_sum),"\nAREA:",
                                  sprintf("%0.2f",area_normalized)),col=vcol,
           y.intersp=1.5,lwd=2,cex = 0.7,pt.cex=0.7)
  }
  return(TRUE)
}

plotCompound <- function(raw,peaktable,name,mz,rt,rt_win=1,mz_win=0.01,addLine=TRUE,prefix="",addLeg=""){
  title <- paste(prefix,name," MZ: ",sprintf("%0.3f",mz),"(+/-",mz_win,") RT: ",sprintf("%0.2f",rt),"(+/-",rt_win,")",sep="")
  vp <- plotPeaks(peaktable,raw,mz_lim=c(mz-mz_win,mz+mz_win),rt_lim=c(rt-rt_win,rt+rt_win),title=title,addLeg=addLeg)
  if(addLine &vp){
    abline(v=rt,lty=2,lwd=1.5)
  }
  return(vp)
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
  return(is_plotted)
}


###Errror are stored in the error_mat matrix

serr <- which(error_mat>5,arr.ind=TRUE)
tsamp_error <- table(serr[,2])
# sel_samp <- as.numeric(names(tsamp_error[order(tsamp_error,decreasing=TRUE)[1:10]]))
tcomp_error <- table(serr[,1])
# sel_comp <- as.numeric(names(tcomp_error[order(tcomp_error,decreasing=TRUE)[1:10]]))

serr <- serr[order(serr[,2]),]

pdf("U:/users/Alexis/presentation/Group_meeting_25_11/errors_integrations_no_leg_thin.pdf")
ounique <- as.numeric(names(tsamp_error)[order(tsamp_error,decreasing = TRUE)])

for(ss in ounique[1:40]){
  sserr <- serr[serr[,2]==ss,,drop=FALSE]
  scomp <- sserr[,1]
  samp_name <- colnames(intMatrix)[ss]
  samp_inj_name <- paste(samp_name,"newinj",sep="_")
  inj_path_peaktable <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/res_thin/ADAP/peaktables/",samp_inj_name,".csv",sep="")
  if(file.exists(inj_path_peaktable)){
    samp_name <- samp_inj_name
  }
  path_peaktable <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/res_thin/ADAP/peaktables/",samp_name,".csv",sep="")
  path_raw <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/mzML/",samp_name,".mzML",sep="")
  xr <- xcmsRaw(path_raw)
  pt <- read.table(path_peaktable,sep=",",header=TRUE)
  for(icp in seq_along(scomp)){
    cp <- scomp[icp]
    found_comp <- plotCompound(raw = xr,peaktable = pt,name = compounds$name[vm2[cp]],
                 mz = compounds$M_H_n[vm2[cp]],rt = compounds$est_RT[vm2[cp]],prefix=paste(samp_name," ",sep=""),rt_win = 0.4,
                 addLeg = FALSE)
    # if(found_comp)
    # legend("topright",legend = c(paste("Quant-Mass:",sprintf("%0.2f",quant_neg[sserr[icp,,drop=FALSE]])),
    #                           paste("LCMS-workflow:",sprintf("%0.2f",intMatrix[sserr[icp,,drop=FALSE]]))))
  }
}
dev.off()


samp_name <- "9_SPF_D_c_12"
path_peaktable <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/res_thin/ADAP/peaktables/",samp_name,".csv",sep="")
path_raw <- paste("U:/users/Alexis/data/data_karin/mml_v2/neg/mzML/",samp_name,".mzML",sep="")
xr <- xcmsRaw(path_raw)
pt <- read.table(path_peaktable,sep=",",header=TRUE)


list_group <- extractUniqueMasses(compounds)




###We extract all the EICs for eahc peak tp pllot them simulatneously eventually.
for(i in which(sapply(list_group,length)==3)){
  
  
  plotCompounds(xr,pt,compounds$name[list_group[[i]]],
                compounds$M_H_n[list_group[[i]]],compounds$est_RT[list_group[[i]]],mz_win=0.01,rt_win = 0.2)
  
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

vvidx <- sample.int(length(path_peaktables),15)


path_raw <- all_samples
cEICs <- mapply(path_peaktables[vvidx],path_raw[vvidx],FUN = exEICs,
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


##pos_pe
list_peaktables <- list.files("U:/users/Alexis/data/data_karin/mml_v2/neg/res/ADAP/peaktables",full.names=TRUE)



pdf("U:/users/Alexis/data/data_karin/mml_v2/tout.pdf")
for(n_compounds in 1:nrow(compounds)){
  plotCompounds(xr,pt,compounds[n_compounds,1],compounds[n_compounds,"M_H_n"],compounds[n_compounds,"est_RT"])
}
dev.off()

####Soe  number of detected features
num_detect_xcms <- apply(dm_xcms[,10:470],1,function(x){sum(!is.na(x))})
num_detect_workflow <- apply(dm[,(10+12):(470+12)],1,function(x){sum(!is.na(x))})

df_num_detection <- data.frame(num_detection=c(num_detect_xcms,num_detect_workflow),
                               software=c(rep("XCMS",length(num_detect_xcms)),rep("workflow",length(num_detect_workflow))))

ggbar <- ggplot(df_num_detection,aes(x=num_detection,fill=software))+geom_histogram(position="dodge",bins=15)+scale_y_log10()+xlab("Number of detection")+
  ylab("log10(count)")+theme(legend.title = element_text(size=16),legend.text=element_text(size=16))+labs(fill="Software")+theme(text = element_text(size=16))+
  guides(fill=FALSE)
png("U:/users/Alexis/presentation/Group_meeting_25_11/barplot_num_detection.png")
plot(ggbar)
dev.off()

###Calculating all CV
sel_qcs_xcms <- which(str_detect(colnames(dm_xcms),fixed("QC")))
sel_qcs_work <- which(str_detect(colnames(dm),fixed("QC")))

cv_xcms <- apply(dm_xcms[,sel_qcs_xcms],1,function(x){
  if(sum(!is.na(x))<2) return(NA)
  sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
})
cv_work <- apply(dm[,sel_qcs_work],1,function(x){
  if(sum(!is.na(x))<2) return(NA)
  sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)
})


fun_count2 <- function(x){
  med <- median(x)
  return(data.frame(label=paste("count:",length(x),"\nmedian:",sprintf("%0.1f",med*100),
                                "\nmean:",sprintf("%0.1f",mean(100*x)),"%",sep=""),y=median(x)+sd(x)*1.3))
}


df_cvs <- data.frame(cv=c(cv_xcms,cv_work),
                     software=c(rep("XCMS",length(cv_xcms)),rep("Workflow",length(cv_work))))
ggp_b <- ggplot(df_cvs,aes(y=cv,x=df_cvs$software,fill=df_cvs$software))+geom_boxplot()+
  theme(text = element_text(size=16))+xlab("Software")+ylab("CV")+labs(fill="Software")+ylim(c(0,2.5))+
  stat_summary(fun.data = fun_count2, geom = "text", fun.y = median,position = position_dodge(width = 0.75),col="black")+
  theme(
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )
png("U:/users/Alexis/presentation/Group_meeting_25_11/barplot_cvs_qcs.png",width=700,height=500)
plot(ggp_b)
dev.off()