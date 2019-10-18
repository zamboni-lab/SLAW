suppressMessages(library(stringr, warn.conflicts = FALSE))
suppressMessages(library(DBI, warn.conflicts = FALSE))
suppressMessages(library(RSQLite, warn.conflicts = FALSE))
suppressMessages(library(CAMERA, warn.conflicts = FALSE))
suppressMessages(library(Matrix, warn.conflicts = FALSE))
suppressMessages(library(MSnbase, warn.conflicts = FALSE))
suppressMessages(library(BiocParallel, warn.conflicts = FALSE))
suppressMessages(library(jsonlite, warn.conflicts = FALSE))
suppressMessages(library(cliqueMS, warn.conflicts = FALSE))
suppressMessages(library(igraph, warn.conflicts = FALSE))
suppressMessages(library(igraph, warn.conflicts = FALSE))
suppressMessages(library(InterpretMSSpectrum, warn.conflicts = FALSE))

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

getIntensityPos <- function(dm){
  which(startsWith(colnames(dm),"intensity"))
}


convertToCliqueMS <- function(dm,
                              path_raw,
                              ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData") {
  suppressMessages(library(MSnbase, warn.conflicts = FALSE))
  suppressMessages(library(xcms, warn.conflicts = FALSE))
  mzdata <- readRDS(file = ref_xcms)

  getIntensityPos <- function(dm){
    which(startsWith("intensity",colnames(dm)))
  }

  pint <- getIntensityPos(dm)
  intensity <- apply(dm[, pint], 1, mean, na.rm = TRUE)

  ###We mdoify all the fileds of the available object
  cnames <- c(
    "mz",
    "mzmin",
    "mzmax",
    "rt",
    "rtmin",
    "rtmax",
    "into",
    "intb",
    "maxo",
    "sn",
    "sample",
    "is_filled"
  )

  mzraw <- readMSData(path_raw,mode="onDisk")
  rtime <- rtime(mzraw)
  tdf <-
    data.frame(
      mz = dm[, "mz"],
      mzmin = dm[, "mz_min"] - 0.0003,
      mzmax = dm[, "mz_max"] + 0.0003,
      rt = 60 * (dm[, "rt"]),
      rtmin = 60 * (dm[, "rt_min"] - dm[,"mean_peakwidth"]/2-0.002),
      rtmax = 60 * (dm[, "rt_max"] + dm[,"mean_peakwidth"]/2+0.002),
      into = intensity,
      intb = intensity,
      maxo = intensity,
      sn = rep(10, nrow(dm)),
      sample =
        rep(1, nrow(dm)),
      is_filled = rep(0, nrow(dm))
    )

  ###We look for the closest retention time in these data
  rttime <- rtime(mzraw)
  blim <- c(-1,rttime[2:length(rttime)]-diff(rttime)/2,100000)
  rtmatch_min <- .bincode(tdf[,"rtmin"],breaks = blim,include.lowest = TRUE)
  rtmatch_min <- ifelse(rtmatch_min>1,rtmatch_min-1,rep(1,length(rtmatch_min)))
  rtmatch_max <- .bincode(tdf[,"rtmax"],breaks = blim,include.lowest = TRUE)
  rtmatch_max <- ifelse(rtmatch_max<length(rttime),rtmatch_max+1,rep(length(rttime),length(rtmatch_max)))
  rtmatch_med <- .bincode(tdf[,"rt"],breaks = blim,include.lowest = TRUE)

  # sel <- apply(rtmatch_min,rtmatch_max)

  tdf[,"rtmin"] <- rttime[rtmatch_min]
  tdf[,"rtmax"] <- rttime[rtmatch_max]
  tdf[,"rt"] <- rttime[rtmatch_med]

  ###only keep feature which have a different rt,rtmin and rtmax
  vdup <- apply(tdf[,c("rtmin","rt","rtmax")],1,anyDuplicated)
  sel <- which(vdup==0)

  tenv <- new.env()
  tenv[["chromPeaks"]] <- as.matrix(tdf[sel,])
  mzdata@msFeatureData@.xData <- tenv
  mzdata@featureData@data <- mzraw@featureData@data
  mzdata@processingData@files <- path_raw
  return(list(mzdata,sel))
}


computeNetworkRawfile <-
  function(dm,
           raw_data,
           ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData") {
    suppressMessages(library(cliqueMS, warn.conflicts = FALSE))
    ####convert a datamatrix to a peaktable
    convertToCliqueMS <- function(dm,
                                  path_raw,
                                  ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData") {
      suppressMessages(library(MSnbase, warn.conflicts = FALSE))
      suppressMessages(library(xcms, warn.conflicts = FALSE))
      mzdata <- readRDS(file = ref_xcms)
      intensity <- apply(dm[, 3:ncol(dm)], 1, mean, na.rm = TRUE)

      ###We change all the ifllled values
      cnames <- c(
        "mz",
        "mzmin",
        "mzmax",
        "rt",
        "rtmin",
        "rtmax",
        "into",
        "intb",
        "maxo",
        "sn",
        "sample",
        "is_filled"
      )

      mzraw <- readMSData(path_raw,mode="onDisk")
      rtime <- rtime(mzraw)
      tdf <-
        data.frame(
          mz = dm[, "mz"],
          mzmin = dm[, "mz"] - 0.0003,
          mzmax = dm[, "mz"] + 0.0003,
          rt = 60 * (dm[, "rt"]),
          rtmin = 60 * (dm[, "rt"] - 0.002),
          rtmax = 60 * (dm[, "rt"] + 0.002),
          into = intensity,
          intb = intensity,
          maxo = intensity,
          sn = rep(10, nrow(dm)),
          sample =
            rep(1, nrow(dm)),
          is_filled = rep(0, nrow(dm))
        )
      ###We look for the closest retention time in these data
      rttime <- rtime(mzraw)
      blim <- c(-1,rttime[2:length(rttime)]-diff(rttime)/2,100000)
      rtmatch_min <- .bincode(tdf[,"rtmin"],breaks = blim,include.lowest = TRUE)
      rtmatch_min <- ifelse(rtmatch_min>1,rtmatch_min-1,rep(1,length(rtmatch_min)))
      rtmatch_max <- .bincode(tdf[,"rtmax"],breaks = blim,include.lowest = TRUE)
      rtmatch_max <- ifelse(rtmatch_max<length(rttime),rtmatch_max+1,rep(length(rttime),length(rtmatch_max)))
      rtmatch_med <- .bincode(tdf[,"rt"],breaks = blim,include.lowest = TRUE)

      # sel <- apply(rtmatch_min,rtmatch_max)

      tdf[,"rtmin"] <- rttime[rtmatch_min]
      tdf[,"rtmax"] <- rttime[rtmatch_max]
      tdf[,"rt"] <- rttime[rtmatch_med]


      ###only keep feature which have a different rt,rtmin and rtmax
      vdup <- apply(tdf[,c("rtmin","rt","rtmax")],1,anyDuplicated)
      sel <- which(vdup==0)

      tenv <- new.env()
      tenv[["chromPeaks"]] <- as.matrix(tdf[sel,])
      mzdata@msFeatureData@.xData <- tenv
      mzdata@featureData@data <- mzraw@featureData@data
      mzdata@processingData@files <- path_raw
      return(list(mzdata,sel))
    }
    ldata <-
      convertToCliqueMS(
        dm,
        path_raw = raw_data,
        ref_xcms = ref_xcms
      )

    mzdata <- ldata[[1]]
    ###We just remove all the EICs value for simplicit


    ##We only extract EICs for feture which seems to fall in the correct boundaries
    sel <- ldata[[2]]
    anclique <- createanClique(mzdata)

    netlist <-
      createNetwork(
        mzdata,
        anclique@peaklist[sel,],
        filter = FALSE,
        mzerror = 1e-5,
        intdiff = 1e4,
        rtdiff = 1e-3
      )$network

    if(is.na(netlist)){
      return(matrix(nrow=0,ncol=3))
    }
    alle <- as_data_frame(netlist, "edges")
    vid <- vertex_attr(netlist, name = "id")
    alle[, 1] <- sel[vid[alle[, 1]]]
    alle[, 2] <- sel[vid[alle[, 2]]]
    alle <- as.matrix(alle)
    return(alle)
  }



createNetworkMultifiles <-
  function(dm,
           raw_files,
           size_batch = 10,
           ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData",
           bpp = NULL) {

    size_batch <- min(size_batch,length(raw_files)-1)
    if (is.null(bpp))
      bpp <- bpparam()

    ###The sparse matrix which will be used ot create the network.
    countMat <- Matrix(0,
                       nrow = nrow(dm),
                       ncol = nrow(dm),
                       sparse = TRUE)
    cosMat <- Matrix(0,
                     nrow = nrow(dm),
                     ncol = nrow(dm),
                     sparse = TRUE)

    ###We build a sequnce with each matrix
    seq_cut <- seq(0, length(raw_files), by = size_batch)
    if (seq_cut[length(seq_cut)] != length(raw_files)) {
      seq_cut <- c(seq_cut, length(raw_files) + 1)
    } else{
      seq_cut[length(seq_cut)] <- length(raw_files) + 1
    }
    message("Processing segment: ")

    for (i in 1:(length(seq_cut) - 1)) {
      message(i, " ",appendLF = FALSE)
      ###We compute the network for the selected files.
      ledges <-
        bplapply(
          as.list(raw_files[seq_cut[i]:(seq_cut[i + 1] - 1)]),
          FUN = computeNetworkRawfile,
          ref_xcms = ref_xcms,
          dm = dm,
          BPPARAM = bpp
        )


      ###We add the edges to the main network
      for (j in seq_along(ledges)) {
        alle <- ledges[[j]]
        if(nrow(alle)==0) next
        # vid <- vertex_attr(ledges[[j]], name = "id")
        # alle[, 1] <- vid[alle[, 1]]
        # alle[, 2] <- vid[alle[, 2]]
        alle <- as.matrix(alle)
        ##We update the cost matrix
        cosMat[alle[, c(1, 2)]] <-
          (cosMat[alle[, 1:2]] * countMat[alle[, 1:2]] + alle[, 3]) /
          (countMat[alle[, 1:2]] + 1)
        cosMat[alle[, c(2, 1)]] <- cosMat[alle[, c(1, 2)]]
        countMat[alle[, c(1, 2)]] <- countMat[alle[, c(1, 2)]] + 1
      }
    }

    message("Done")

    ##We get the first useless elemnts
    sel_val <- rowSums(cosMat)!=0
    # sel_val <- apply(cosMat, 1, function(x) {
    #   any(x != 0)
    # })
    cosMat <- cosMat[sel_val, sel_val]

    gadj <-
      graph_from_adjacency_matrix(adjmatrix = cosMat,
                                  mode = "undirected",
                                  weighted = TRUE)
    gadj <-
      set_vertex_attr(gadj, name = "id", value = which(sel_val))
    ldata <-
      convertToCliqueMS(
        dm,
        path_raw = raw_files[[1]],
        ref_xcms = ref_xcms
      )

    mzdata <- ldata[[1]]

    intensity <- apply(dm,1,function(x){mean(x[3:length(x)],na.rm=TRUE)})
    anclique <- createanClique(mzdata)
    anclique$peaklist <- data.frame(
      mz = dm[, "mz"],
      mzmin = dm[, "mz"] - 0.0003,
      mzmax = dm[, "mz"] + 0.0003,
      rt = 60 * (dm[, "rt"]),
      rtmin = 60 * (dm[, "rt"] - 0.002),
      rtmax = 60 * (dm[, "rt"] + 0.002),
      into = intensity,
      intb = intensity,
      maxo = intensity,
      sn = rep(10, nrow(dm)),
      sample =
        rep(1, nrow(dm)),
      is_filled = rep(0, nrow(dm))
    )

    # browser()



    anclique$network <- gadj
    return(anclique)
  }


###Annotations of cliques
annotateCliques <- function(cliques, adducts, main_adducts,
                            ionization_mode,val_int, bpp = NULL) {
  if (is.null(bpp))
    bpp <- bpparam()

  # vfeat <-
  #   bplapply(
  #     cliques$cliques,
  #     FUN = annotateCliqueInterpretMSspectrum,
  #     BPPARAM = bpp,
  #     dm = cliques$peaklist,adducts=adducts,main_adducts=main_adducts,
  #     ionization_mode=ionization_mode
  #   )

  vfeat <-
    bplapply(
      cliques$cliques,FUN = annotateCliqueInterpretMSspectrum,
      dm = cliques$peaklist,adducts=adducts,main_adducts=main_adducts,
      ionization_mode=ionization_mode,val_int = val_int,BPPARAM = bpp
    )

  # vfeat <- lapply(cliques$cliques[1:maxCliques],annotateCliqueIntepretMSspectrum,dm=cliques$peaklist)
  totfeat <- 1
  for (i in seq_along(vfeat)) {
    seq_idx <-
      rep(totfeat:(totfeat + length(vfeat[[i]]) - 1), times = sapply(vfeat[[i]], nrow))
    vfeat[[i]] <- do.call(rbind, vfeat[[i]])
    vfeat[[i]] <- cbind(vfeat[[i]], seq_idx)
    totfeat <- totfeat + nrow(vfeat[[i]])

  }
  ###We merge all the data matrix
  vfeat <- do.call(rbind, vfeat)
  ###We add an idx group to each feature we add the rest of the data
  cnames <-
    c(
      "index",
      "mz",
      "intb",
      "isogr",
      "iso",
      "charge",
      "adduct",
      "ppm",
      "label",
      "group_label"
    )
  colnames(vfeat) <- cnames
  return(vfeat)
}


annotateCliqueInterpretMSspectrum <-
  function(clique,adducts,main_adducts, ionization_mode,
           dm,val_int, ppm = 10, dmz = 0.005) {
    library(InterpretMSSpectrum)

    if (length(clique) == 1) {
      return(list(
        data.frame(
          index = clique[1],
          mz = dm[clique[1], 1],
          intb = val_int[clique[1]],
          isogr = NA,
          iso = NA,
          charge = NA,
          adduct = "[M+H]+",
          ppm = NA,
          label = "[M+H]+"
        )
      ))
    }
    sel_clust_idx <- clique
    sel_idx <- seq_along(sel_clust_idx)

    current_val <- cbind(dm[sel_clust_idx,1],val_int[sel_clust_idx])
    colnames(current_val) <- c("mz","int")
    # message("dd_", paste(dim(current_val),collapse = "=="))

    all_features <- vector(mode = "list", length = 10)
    num_feat <- 1

    while (TRUE) {
      ###We find all the annotated adducts.
      if (length(sel_idx) == 1) {
        all_features[[num_feat]] <- data.frame(
          index = sel_clust_idx[sel_idx],
          mz = dm[clique[sel_idx], 1],
          intb = val_int[clique[sel_idx]],
          isogr = NA,
          iso = NA,
          charge = NA,
          adduct = "[M+H]+",
          ppm = NA,
          label = "[M+H]+"
        )
        num_feat <- num_feat + 1
        break
      }
      # cat("INMAIN")
      # print(current_val[sel_idx,])
      # cat("dim",dim(current_val[sel_idx,]))
      annots <-
        findMAIN(
          current_val[sel_idx,,drop=FALSE],ionmode = ionization_mode,
          rules=adducts,adducthyp = main_adducts,
          mzabs = dmz,ppm = ppm,mainpkthr = 0.2
        )[[1]]
      # cat("OUTMAIN")
      sel_adducts <- which(!is.na(annots[, "adduct"]))
      if (length(sel_adducts) <= 0)
        break

      sel_iso <- annots$isogr[sel_adducts]
      sel_iso <- sel_iso[!is.na(sel_iso)]

      sel_iso <- which(annots$isogr %in% sel_iso)
      sel_feat <- union(sel_adducts, sel_iso)
      # if(14434 %in% sel_clust_idx[sel_idx[sel_feat]]) browser()

      all_features[[num_feat]] <-
        cbind(sel_clust_idx[sel_idx[sel_feat]], annots[sel_feat,])
      colnames(all_features[[num_feat]]) <-
        c("index",
          "mz",
          "intb",
          "isogr",
          "iso",
          "charge",
          "adduct",
          "ppm",
          "label")
      num_feat <- num_feat + 1
      ##We remove the selected features.
      sel_idx <- sel_idx[-sel_feat]

    }
    all_features <- all_features[1:(num_feat - 1)]
    return(all_features)
  }

convertFeatures <- function(resAnnot){
  ####We change all the data
  unique_groups <- unique(resAnnot$group_label)
  resList <- vector(mode="list",length = length(unique_groups))
  current_dec <- 0
  message("Features conversion: ")
  for(ig in seq_along(resList)){
    if(((ig*10) %/% length(resList))!=current_dec){
      message(current_dec*10," ",appendLF = FALSE)
      current_dec <- (ig*10) %/% length(resList)
    }
    xannot <- resAnnot[resAnnot$group_label==unique_groups[ig],]
    vannot <- unique(xannot[,"isogr"])
    vannot <- vannot[!is.na(vannot)]
    p_unique <- NULL
    if(length(vannot)!=0){
      p_unique <- match(vannot,xannot[,"isogr"])
    }else{
      p_unique <- match(vannot,xannot[,"isogr"])
    }

    sel_adducts <- xannot[,"label"]
    vvv <- which(!is.na(xannot[,"adduct"]))
    pmax <- vvv[which.max(xannot[vvv,"intb"])]
    vs <- InterpretMSSpectrum:::getRuleFromIonSymbol(xannot[pmax,"label"])
    neutral_mass <- as.numeric((vs[2]*xannot[pmax,"mz"]-vs[4])/vs[3])
    tempdf <- data.frame(index=xannot[,"index"],neutral_mass=rep(neutral_mass,nrow(xannot)),label=xannot[,"label"],
                         charge=xannot[,"charge"],group_label=xannot[,"group_label"],
                         adduct=xannot[,"label"],ref_feature=rep(xannot[pmax,"index"],nrow(xannot)),
                         stringsAsFactors = FALSE)
    if(length(vannot)==0){
      resList[[ig]] <- tempdf
      next
    }
    for(p in seq_along(p_unique)){
      pp <- p_unique[p]
      pisos <- which(xannot[,"isogr"]==xannot[pp,"isogr"])
      ref_label <- xannot[pp,"label"]
      if(length(pisos)==1){
        tempdf[pisos,"adduct"] <- ref_label
        next
      }
      if(is.na(ref_label)) ref_label <- "[M+H]+"

      ###We update the label of all the peaks
      labelC13 <- paste(ref_label,c("",paste("+",1:(length(pisos)-1),sep="")))
      tempdf[pisos,"adduct"] <- labelC13
    }
    ###We reorder the data.frame
    tempdf <- tempdf[order(tempdf$index),]
    resList[[ig]] <- tempdf
  }
  return(resList)
}


###The full data matrix annotated.
buildDataMatrixFull <- function(dm,annot){
  vannot <- sapply(annot,function(an,dm,cnames){
    sub_dm <- cbind(dm[an[,"index"],c(1,2)],an[,"adduct"],an[,"group_label"],
          an[,"neutral_mass"],dm[an[,"index"],3:ncol(dm)])

    colnames(sub_dm) <- c("mz","rt","annotation","group","neutral_mass",cnames[3:ncol(dm)])
    return(sub_dm)
  },dm=dm,cnames=colnames(dm),
  simplify=FALSE)

  vannot <- do.call(rbind,vannot)
  return(vannot)
}


###A single line by data matrix
buildDataMatrixSimplified <- function(dm,annot){
  vannot <- sapply(annot,function(an,dm,cnames){
    refv <- an[1,"ref_feature"]
    main_adduct <- an[match(refv,an[,1]),"adduct"]
    all_adducts <- paste(an[,"adduct"],collapse = "|")
    neutral_mass <- an[1,"neutral_mass"]
    tdf <- data.frame(mz=dm[refv,1],rt=dm[refv,2],main_adduct=main_adduct,all_adducts=all_adducts,
                      neutral_mass=neutral_mass)
    tdf <- cbind(tdf,dm[refv,3:ncol(dm)])

    colnames(tdf) <- c("mz","rt","main_peak","annotations","neutral_mass",cnames[3:ncol(dm)])
    # browser()
    return(tdf)
  },dm=dm,cnames=colnames(dm),
  simplify=FALSE)

  vannot <- do.call(rbind,vannot)
  return(vannot)
}


groupFeatures <- function(dm,val_int, raw_files, adducts,main_adducts,ionization_mode,
                          ppm = 10, dmz = 0.005, size_batch = 10,
                          ref_xcms="X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData", bbp=NULL){
  if(is.null(bpp)) bpp <- bpparam()

  anclique <- createNetworkMultifiles(dm,
             raw_files,
             size_batch = size_batch,
             ref_xcms = ref_xcms,
             bpp = bpp)
  ###we compute the cliques
  anclique <- computeCliques(anclique, 1e-5, TRUE)
  pint <- getIntensityPos(dm)[1]

  res_df <- annotateCliques(anclique,adducts,main_adducts,ionization_mode, val_int, bpp = bpp)
  ####We add the annotation for all the peaks
  annot <- convertFeatures(res_df)
  dm_full <- buildDataMatrixFull(dm,annot)
  dm_simple <- buildDataMatrixSimplified(dm,annot)
  return(list(dm_simple,dm_full))
}






###Examples of commandline

###Rscript annotation_mixed_method.R X:/Documents/dev/LCMSfeatureGrouping/inst/datamatrices_ecoli_philipp/pos_ecoli_splash_mix_group_2.csv

###PATH_DATAMATRIX

args <- commandArgs(trailingOnly = TRUE)

###for testing purpose

# args <- c("X:/Documents/dev/LCMSfeatureGrouping/inst/datamatrices_ecoli_philipp/pos_ecoli_splash_mix_group_2.csv",
#           "U:/users/Alexis/data/philipp_ecoli_splash_mix/res_pos/pos_splash_mix_ecoli_local.db","test_mat_simple.csv",
#           "test_mat_fule.csv","3","X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData",
#           "C:/Users/dalexis/Documents/python/pylcmsprocessing/data/adducts_pos.txt",
# "C:/Users/dalexis/Documents/python/pylcmsprocessing/data/adducts_main_pos.txt","positive","110",
# "0.005","2")
###PATH_DATAMATRIX


PATH_DATAMATRIX <- args[1]
PATH_DB <- args[2]
PATH_OUTPUT_FULL <- args[3]
PATH_OUTPUT_SIMPLE <- args[4]
NUM_CORES <- as.numeric(args[5])
PATH_MODEL <- args[6]
PATH_ADDUCTS <- args[7]
PATH_MAIN_ADDUCTS <- args[8]
POLARITY <- args[9]
PPM <-  as.numeric(args[10])
DMZ <-  as.numeric(args[11])
FILTER_NUMS <- max(1,as.numeric(args[13]))

###We peak the FILE_USED most intense files.
FILES_USED <- as.numeric(args[12])

###reading data matrices
dm <- read.table(PATH_DATAMATRIX,header=TRUE,sep=",")

posIntensities <- getIntensityPos(dm)

# message(paste(posIntensities,"|"))

# message(paste(posIntensities,collapse = "|"))


num_detect <- apply(dm[,posIntensities,drop=FALSE],1,function(x){sum(!is.na(x))})
message(paste(as.numeric(table(num_detect)),collapse = ".."))
num_detect <- num_detect>=FILTER_NUMS

message("Retained ",sum(num_detect),"files on ",nrow(dm))
dm <- dm[num_detect,,drop=FALSE]



###Reading the raw files
dbb <- dbConnect(RSQLite:::SQLite(),PATH_DB)
raw_files <- dbGetQuery(dbb,"SELECT path FROM samples")[,1]
dbDisconnect(dbb)

####Selecting the msot intense files
val_int <- apply(dm[,posIntensities],2,sum,na.rm=TRUE)
sel_files <- order(val_int,decreasing = TRUE)[1:min(FILES_USED,length(val_int))]
raw_files <- raw_files[sel_files]


val_int_var <- apply(dm[,posIntensities],1,mean,na.rm=TRUE)
###Setting up the parallel processing
bpp <- NULL
if(get_os()=="win"){
  bpp <- SnowParam(workers = NUM_CORES)
}else{
  bpp <- MulticoreParam(workers = NUM_CORES)
}
# bpp <- SerialParam()

fadd <- file(PATH_ADDUCTS,"r")
adducts <- readLines(fadd)
close(fadd)

###
fadd <- file(PATH_MAIN_ADDUCTS,"r")
main_adducts <- readLines(fadd)
close(fadd)


# annotated_tables <- groupFeatures(dm[sort(sample.int(nrow(dm),size=5000)),], raw_files[1:2], adducts,main_adducts,ionization_mode="positive",
#                                   mzwin=MZWIN,rtwin=RTWIN,ppm = PPM, dmz = DMZ, size_batch = NUM_CORES,
#                                   ref_xcms=PATH_MODEL, bbp=bpp)

annotated_tables <- groupFeatures(dm,val_int_var, raw_files, adducts,main_adducts,ionization_mode=POLARITY,
                          ppm = PPM, dmz = DMZ, size_batch = NUM_CORES,
                          ref_xcms=PATH_MODEL, bbp=bpp)

write.table(annotated_tables[[1]],file = PATH_OUTPUT_SIMPLE,sep=";",row.names = FALSE)

write.table(annotated_tables[[2]],file = PATH_OUTPUT_FULL,sep=";",row.names = FALSE)

###for testing purpose
