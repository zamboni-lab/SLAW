###Loading the DB connection of samples
suppressWarnings(suppressMessages(library(RSQLite, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(DBI, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(stringr, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(MSnbase, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(igraph, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(rtree, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))

sink(file=stdout())
###Merging of spectra using the msClust algorithm
mergeSpectra <- function(x,specs,tab_summary,cos_thresh=0.7,round=10){

    ###We always bin the spectra
    binv <- seq(0.5,2000.5,by=0.1)

    approxCos <- function(s1,s2,npeaks=2){
        v1 <- .bincode(s1[,1],breaks = binv)
        v2 <- .bincode(s2[,1],breaks = binv)

        f1 <- s1[,1]**2*s1[,2]**0.5
        f2 <- s2[,1]**2*s2[,2]**0.5

        ###We compute the intersection
        vmm <- match(v1,v2)

        common_s1 <- which(!is.na(vmm))
        common_s2 <- vmm[common_s1]

        if(length(common_s1)<npeaks) return(0)

        only_s1 <- which(is.na(vmm))
        only_s2 <- setdiff(seq_along(f2),common_s2)

        ###We penalize it by the intensity of the non matche dpeaks
        cos_term <- sum(f1[common_s1]*f2[common_s2])/
        (sqrt(sum(f1[common_s1]**2))*sqrt(sum(f2[common_s2]**2)))

        if(length(only_s2)==0 & length(only_s1)==0) return(cos_term)

        0.5*cos_term +0.5*(sum(f1[common_s1]+f2[common_s2]))/(sum(f1)+sum(f2))
    }

    fastCluster <- function(lspec,tau_min=0.5,nrounds=20){
        if(length(lspec)==1) return(list(list(1),list(lspec),1))
        change <- TRUE
        tau <- 0.999
        delta <- (1-tau_min)/nrounds
        vclust <- seq_along(lspec)
        to_consider <- rep(TRUE,length(lspec))

        ###Actual informatiosn about clusters.
        clust_idx <- as.list(seq_along(lspec))
        clust_size <- rep(1,length(lspec))
        representative_spec <- lspec
        representative_idx <- seq_along(lspec)
        ncluster <- length(lspec)
        for(i in 1:nrounds){
            if(!change) break
            change <- FALSE
            sel_idx <- which(to_consider)
            if(length(sel_idx)<=1) break
            ssel_idx <- seq_along(sel_idx)
            ssel_idx <- 2:length(sel_idx)
            for(ic in ssel_idx){
                min_vec <- sel_idx
                sic <- sel_idx[ic]
                if(is.na(sic)) next
                sco <- to_consider[1:ic]
                if(sum(sco)==0) next
                sub_idx <- (1:(ic-1))[to_consider[1:(ic-1)]]
                for(iic in seq_along(sub_idx)){
                    siic <- sub_idx[iic]
                    if(is.na(siic)) next
                    if(approxCos(representative_spec[[sic]],representative_spec[[siic]])>=tau){
                        ###We merge the cluster
                        clust_idx[[siic]] <- c(clust_idx[[siic]],clust_idx[[sic]])
                        ###We keep the spectrum with the highest number of peaks.
                        if(nrow(representative_spec[[sic]])>nrow(representative_spec[[siic]])){
                            representative_spec[[siic]] <- representative_spec[[sic]]
                            representative_idx[siic] <- representative_idx[sic]
                        }
                        ##We update the cluster in all cases
                        clust_size[siic] <- clust_size[siic]+clust_size[sic]
                        sel_idx[sic] <- NA
                        to_consider[sic] <- FALSE
                        change <- TRUE
                        break
                    }
                }
            }
            tau <- tau-tau_min
        }
        ###Returning the cluster and the idx
        return(list(clust_idx[to_consider],representative_spec[to_consider],
        representative_idx[to_consider]))
    }


    ###extracting the data.frame
    specs <- apply(tab_summary[x,,drop=FALSE],1,function(x,ref){ref[[x[4]]][[x[1]]]},ref=specs)
    specs <- sapply(specs,as.data.frame,simplify=FALSE)
    fc <- fastCluster(specs,nrounds=round,tau_min = cos_thresh)

    ##We keep the biggest components, or the one with the most peaks
    size_comp <- sapply(fc[[1]],length)
    ppmax <- which.max(size_comp)
    ####We then return the informations, number of spectra avraged, number of spectra discarded
    return(list(spec=fc[[2]][[ppmax]],idx=x[fc[[3]][ppmax]],num_total=length(x),num_fused=size_comp[ppmax]))

}


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


###Modified Mgf Writing function
writeMgfDataFileToConnection <- function (splist, con, TITLE = NULL,
addFields = NULL)
{
    if (length(addFields)) {
        if (length(dim(addFields)) != 2)
        stop("'addFields' has to be a matrix or data.frame.")
        if (!is.matrix(addFields))
        addFields <- do.call(cbind, lapply(addFields, as.character))
        if (is.null(colnames(addFields)))
        stop("Column names required on 'addFields'.")
        if (nrow(addFields) != length(splist))
        stop("nrow of 'addFields' has to match length of 'splist'")
    }
    for (i in seq(along = splist)) {
        MSnbase:::writeMgfContent(splist[[i]], TITLE = TITLE, con = con,
        addFields = addFields[i, ])
    }
}



args <- commandArgs(trailingOnly = TRUE)
PATH_DB <- args[1]
MZ_TOL <- as.numeric(args[2])
RT_TOL <- as.numeric(args[3])
NUM_CORES <- as.numeric(args[4])
PATH_MGF <- args[5]
TEMP_LOCATION <- args[6]
TEMP_FILE_SWITCHING <- args[7]


##Can be changed
BY_BATCH <- 7000
bpp <- NULL

if (get_os() == "win") {
    bpp <- SnowParam(workers = NUM_CORES,progressbar=FALSE)
} else{
    bpp <- MulticoreParam(workers = min(NUM_CORES, 4),progressbar=FALSE)
}

##We get the datamatrxi path
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
PATH_DATAMATRIX <- dbGetQuery(dbb, "SELECT peaktable FROM peakpicking")[, 1]
dbDisconnect(dbb)

##We get all the MS2 paths.
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
PATH_MSMS <- dbGetQuery(dbb, "SELECT output_ms2 FROM processing")[, 1]
dbDisconnect(dbb)

# Some file can be missed if a file is corrupted, so we remove non existing file.
PATH_MSMS <- PATH_MSMS[file.exists(PATH_MSMS)]
### Checking if an actual file exist
if(length(PATH_MSMS)==0) stop("No MS2 spectra to Merge.")

###We read the data matrix.
dmm <- fread(PATH_DATAMATRIX,sep="\t",header=TRUE)
dmm <- dmm[,c("mz","rt","num_detection"),drop=FALSE]
dmm <- as.data.frame(dmm)
num_line <- nrow(dmm)

##Rttol is 0.5 % of your datasets.
if(is.na(RT_TOL)){
    RT_TOL <- diff(range(dmm[,"rt"]))/200
}

##We create an Rtree to map the features to it.
vmat <- as.matrix(dmm[,c("rt","mz")])
vmat[,1] <- vmat[,1]/RT_TOL
vmat[,2] <- vmat[,2]/MZ_TOL
rtt <- RTree(vmat)

####We read the mgf data from the disk.
vmgf <- bplapply(PATH_MSMS,readMgfData,BPPARAM = bpp)

#### Building a summary table
tab_summary <- sapply(vmgf,function(x){
  header(x)[,c(2,3)]
  },simplify=FALSE)

###We extract the collision energy 
collision <- NULL
if("COLLISIONENERGY" %in% fvarLabels(vmgf[[1]])){
  collision <- sapply(vmgf,function(x){as.numeric(as.character(fData(x)[,"COLLISIONENERGY"]))},simplify = FALSE)
}else{
  ###In this case we don t know the collision energy, we set it to 0
  collision <- sapply(vmgf,function(x){rep(0,length(x))},simplify = FALSE)
}

collision <- unlist(collision)

tab_summary <- mapply(tab_summary,seq_along(vmgf),FUN=function(x,y){
    cnames <- colnames(x)
    x <- cbind(1:nrow(x),x,rep(y,nrow(x)))
    colnames(x) <- c("idx",cnames,"file")
    return(x)
},SIMPLIFY = FALSE)

###We look for the nearest neighbour for each peak using the eclidien distnce weighted by the tolerance
tab_summary <- do.call(rbind,tab_summary)
tab_summary[,"precursor.mz"] <- tab_summary[,"precursor.mz"]/MZ_TOL
tab_summary[,"retention.time"] <- tab_summary[,"retention.time"]/(60*RT_TOL)


####We find the nearest neighbour for eveyr spectrum.
vclose <- rtree::knn.RTree(rtt,as.matrix(tab_summary[,c("retention.time","precursor.mz")]),as.integer(1))
vbox <- rtree::withinBox.RTree(rtt,as.matrix(tab_summary[,c("retention.time","precursor.mz")]),dx=1,
dy=1)

##We calculate the intersection which give the real matched feature every time.
def_val <- mapply(vclose,vbox,FUN=function(x,y){
    intersect(x,y+1)
},SIMPLIFY = FALSE)

##We only keep the feature with a matching MS-MS spectra.
is_fragmented <- which(sapply(def_val,function(x){length(x)>0}))
if(length(is_fragmented)==0) stop("No MS-MS spectra found with matching LC-MS feature.")
message("Found ",length(is_fragmented)," features with associated MS-MS spectra")

###We aggregate the  data by nearest neighbours and energy
vfactor <- apply(cbind(unlist(def_val[is_fragmented]),collision[is_fragmented]),1,paste,collapse="_")
listidx <- by(is_fragmented,INDICES = vfactor,FUN=function(x){x})

####We now process the data by batch to avoid any overhead
fcc <- bplapply(listidx,FUN = mergeSpectra,specs=vmgf,tab_summary=tab_summary,BPPARAM=bpp)

###We extract all the necessary spectra
dm_idx <- str_split(names(listidx),fixed("_"),simplify = TRUE)
dm_idx <- apply(dm_idx,2,as.numeric)
num_fused <- sapply(fcc,function(x){x[[4]]})
spec_idx <- sapply(fcc,"[[",i=2)

###We build a list for all the spectra.
consensus_specs <- apply(tab_summary[spec_idx,,drop=FALSE],1,function(x,ref){ref[[x[4]]][[x[1]]]},ref=vmgf)


###WE make a table of fileds to add
supp_infos <- data.frame(SCANS=1:nrow(dm_idx),FEATURE=dm_idx[,1],ENERGY=dm_idx[,2],NUM_CLUSTERED=sapply(fcc,function(x){x[[4]]}))

###We store the feature information into a file.
##We always write the spectra
tcon <- file(description = PATH_MGF, open = "w")
writeMgfDataFileToConnection(consensus_specs, con = tcon,
addFields = supp_infos)
close.connection(tcon)

##We now add two oclumns ot the data matrix add two columns to the data matrix to avoid later complication.
###We rewrite the data matrix by batch

###Reading the column name first.
ocnames <- as.character(fread(PATH_DATAMATRIX,sep = "\t",nrows=1,header=FALSE)[1,])

###We detect the position of the first qaunt_columns
quant_prefix <- paste(str_split(ocnames[length(ocnames)],fixed("_"))[[1]][1],"_",sep="")
to_cut <- which(startsWith(ocnames,quant_prefix))[1]
df_meta <- data.frame(feature=dm_idx[,1],energy=dm_idx[,2],index=1:nrow(dm_idx))
id_ener_summary <- by(df_meta,INDICES =df_meta$feature,FUN=function(x){
  paste(x[["index"]],paste("(e",x[["energy"]],")",sep=""),sep="_",collapse = "|")
})
num_fused_all <- tapply(num_fused,INDEX = df_meta$feature,FUN = sum)
pos_dm <- as.integer(names(id_ener_summary))
o_dm_idx <- order(pos_dm,decreasing=FALSE)
first_spec <- last_spec <- 0
seq_cut <- seq(2,nrow(dmm),by=BY_BATCH)
if(seq_cut[length(seq_cut)]!=nrow(dmm)){
    seq_cut[length(seq_cut)+1] <- nrow(dmm)
}

###We then write the file to the data by batch
for(i in 1:(length(seq_cut)-1)){
    firstLine <- seq_cut[i]
    lastLine <- seq_cut[i+1]-1

    ###We recover the spec to write
    first_spec <- last_spec+1
    last_spec <- first_spec

    while((pos_dm[o_dm_idx[last_spec]]<(lastLine-1))&
    (last_spec<=length(o_dm_idx))){
        last_spec <- last_spec+1
    }
    last_spec <- last_spec-1


    sub_dm <- fread(PATH_DATAMATRIX,sep = "\t",skip=firstLine-1,nrows=lastLine-firstLine+1,header=TRUE)
    colnames(sub_dm) <- ocnames

    ###We write the elemnts in the new column
    seq_ms2_idx <- rep(NA,nrow(sub_dm))
    seq_num_ms2 <- rep(NA,nrow(sub_dm))

    ###Adding the supplementary MS2 informations to the table.
    if(last_spec>first_spec){
        ###We always verify that the msms spectra are non negative
        sel_pos <- pos_dm[o_dm_idx[first_spec:last_spec]]
        sel_ppos <- sel_pos>(firstLine-1)
        if(sum(sel_ppos)>=1){
          seq_ms2_idx[pos_dm[o_dm_idx[first_spec:last_spec]][sel_ppos]-firstLine+1] <- id_ener_summary[o_dm_idx[first_spec:last_spec][sel_ppos]]
          seq_num_ms2[pos_dm[o_dm_idx[first_spec:last_spec]][sel_ppos]-firstLine+1] <- num_fused_all[o_dm_idx[first_spec:last_spec][sel_ppos]]
        }
    }

    ###WE add it to the data table eventually
    cnames <- colnames(sub_dm)
    cnames <- c(cnames[1:(to_cut-1)],"ms2_id","num_clustered_ms2",cnames[to_cut:ncol(sub_dm)])
    sub_dm <- cbind(sub_dm[,1:(to_cut-1),drop=FALSE],seq_ms2_idx,seq_num_ms2,
    sub_dm[,to_cut:ncol(sub_dm),drop=FALSE])
    colnames(sub_dm) <- cnames
    if(i==1){
        fwrite(sub_dm,file=TEMP_LOCATION,sep = "\t",row.names = FALSE,col.names = TRUE)
    }else{
        fwrite(sub_dm,file=TEMP_LOCATION,sep = "\t",append=TRUE,row.names = FALSE,col.names = FALSE)
    }
}

###We then rename and remove the file.
a <- file.rename(PATH_DATAMATRIX,TEMP_FILE_SWITCHING)
a <- file.rename(TEMP_LOCATION,PATH_DATAMATRIX)
a <- file.remove(TEMP_FILE_SWITCHING)
