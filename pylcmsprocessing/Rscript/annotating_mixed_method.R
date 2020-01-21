suppressWarnings(suppressMessages(library(stringr, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(DBI, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(RSQLite, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(CAMERA, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(Matrix, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(MSnbase, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(jsonlite, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(cliqueMS, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(igraph, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(igraph, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(Rcpp, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(InterpretMSSpectrum, warn.conflicts = FALSE)))

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

getIntensityPos <- function(dm) {
  which(startsWith(colnames(dm), "int"))
}

convertToCliqueMS <- function(dm,
                              path_raw,
                              mzraw,
                              idx,
                              ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData",
                              maxPeaks = 5000) {
  suppressMessages(library(MSnbase, warn.conflicts = FALSE))
  suppressMessages(library(xcms, warn.conflicts = FALSE))
  mzdata <- readRDS(file = ref_xcms)
  sel_idx <- NULL
  if (missing(idx)) {
    sel_idx <- 1:nrow(dm)
  } else{
    sel_idx <- which(!is.na(dm[, idx]))
  }

  if (length(sel_idx) > maxPeaks &
      (!missing(idx)))
    sel_idx <- sort(sample(sel_idx, maxPeaks))

  getIntensityPos <- function(dm) {
    which(startsWith("int", colnames(dm)))
  }
  # pint <- getIntensityPos(dm)
  intensity <- NULL
  if (missing(idx)) {
    intensity <-
      apply(dm[, getIntensityPos(dm), drop = FALSE], 1, mean, na.rm = TRUE)
  } else{
    intensity <- dm[sel_idx, idx]
  }


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
  rtime <- rtime(mzraw)
  tdf <-
    data.frame(
      mz = dm[sel_idx, "mz"],
      mzmin = dm[sel_idx, "mz"] - (dm[sel_idx, "max_mz"] + dm[sel_idx, "min_mz"]) *
        1.5 ,
      mzmax = dm[sel_idx, "mz"] + (dm[sel_idx, "max_mz"] - dm[sel_idx, "min_mz"]) *
        1.5,
      rt = 60 * (dm[sel_idx, "rt"]),
      rtmin = 60 * (dm[sel_idx, "rt"] - dm[sel_idx, "mean_peakwidth"] -
                      0.002),
      rtmax = 60 * (dm[sel_idx, "rt"] + dm[sel_idx, "mean_peakwidth"] +
                      0.002),
      into = intensity,
      intb = intensity,
      maxo = intensity,
      sn = rep(10, length(sel_idx)),
      sample =
        rep(1, length(sel_idx)),
      is_filled = rep(0, length(sel_idx))
    )

  ###We look for the closest retention time in these data
  rttime <- rtime(mzraw)
  blim <- c(-1000, rttime[2:length(rttime)] - diff(rttime) / 2, 100000)
  rtmatch_min <-
    .bincode(tdf[, "rtmin"], breaks = blim, include.lowest = TRUE)
  rtmatch_min <-
    ifelse(rtmatch_min > 1, rtmatch_min - 1, rep(1, length(rtmatch_min)))
  rtmatch_max <-
    .bincode(tdf[, "rtmax"], breaks = blim, include.lowest = TRUE)
  rtmatch_max <-
    ifelse(rtmatch_max < length(rttime),
           rtmatch_max + 1,
           rep(length(rttime), length(rtmatch_max)))
  rtmatch_med <-
    .bincode(tdf[, "rt"], breaks = blim, include.lowest = TRUE)

  # sel <- apply(rtmatch_min,rtmatch_max)

  tdf[, "rtmin"] <- rttime[rtmatch_min]
  tdf[, "rtmax"] <- rttime[rtmatch_max]
  tdf[, "rt"] <- rttime[rtmatch_med]

  ###only keep feature which have a different rt,rtmin and rtmax
  vdup <- apply(tdf[, c("rtmin", "rt", "rtmax")], 1, anyDuplicated)
  sel <- which(vdup == 0)

  tenv <- new.env()
  tenv[["chromPeaks"]] <- as.matrix(tdf[sel, ])
  mzdata@msFeatureData@.xData <- tenv
  mzdata@featureData@data <- mzraw@featureData@data
  mzdata@processingData@files <- path_raw
  return(list(mzdata, sel, sel_idx))
}


computeNetworkRawfile <-
  function(dm,
           idx,
           raw_data,
           mzraw,
           ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData",
           cosFilter = 0.3,
           maxPeaks = 5000) {
    suppressMessages(library(cliqueMS, warn.conflicts = FALSE))
    suppressMessages(library(igraph, warn.conflicts = FALSE))
    ####convert a datamatrix to a peaktable
    convertToCliqueMS <- function(dm,
                                  idx,
                                  path_raw,
                                  mzraw,
                                  ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData") {
      suppressMessages(library(MSnbase, warn.conflicts = FALSE))
      suppressMessages(library(xcms, warn.conflicts = FALSE))
      sel_idx <- which(!is.na(dm[,idx]))
      mzdata <- readRDS(file = ref_xcms)
      sel_idx <- NULL
      if (missing(idx)) {
        sel_idx <- 1:nrow(dm)
      } else{
        sel_idx <- which(!is.na(dm[, idx]))
      }
      # if needed add the maxPeks argument ot the function
      #if(length(sel_idx)>maxPeaks&(!missing(idx))) sel_idx <- sort(sample(sel_idx,maxPeaks))

      getIntensityPos <- function(dm) {
        which(startsWith("int", colnames(dm)))
      }
      # pint <- getIntensityPos(dm)
      intensity <- NULL
      if (missing(idx)) {
        intensity <-
          apply(dm[, getIntensityPos(dm), drop = FALSE], 1, mean, na.rm = TRUE)
      } else{
        intensity <- dm[sel_idx, idx]
      }


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

      ###We create a fake xcmsSet object.
      # mzraw <- readMSData(path_raw, mode = "onDisk")
      rtime <- rtime(mzraw)
      tdf <-
        data.frame(
          mz = dm[sel_idx, "mz"],
          mzmin = dm[sel_idx, "mz"] - (dm[sel_idx, "max_mz"] + dm[sel_idx, "min_mz"]) *
            1.5 ,
          mzmax = dm[sel_idx, "mz"] + (dm[sel_idx, "max_mz"] - dm[sel_idx, "min_mz"]) *
            1.5,
          rt = 60 * (dm[sel_idx, "rt"]),
          rtmin = 60 * (dm[sel_idx, "rt"] - dm[sel_idx, "mean_peakwidth"] -
                          0.002),
          rtmax = 60 * (dm[sel_idx, "rt"] + dm[sel_idx, "mean_peakwidth"] +
                          0.002),
          into = intensity,
          intb = intensity,
          maxo = intensity,
          sn = rep(10, length(sel_idx)),
          sample =
            rep(1, length(sel_idx)),
          is_filled = rep(0, length(sel_idx))
        )

      ###We look for the closest retention time in these data
      rttime <- rtime(mzraw)
      blim <-
        c(-1000, rttime[2:length(rttime)] - diff(rttime) / 2, 100000)
      rtmatch_min <-
        .bincode(tdf[, "rtmin"], breaks = blim, include.lowest = TRUE)
      rtmatch_min <-
        ifelse(rtmatch_min > 1, rtmatch_min - 1, rep(1, length(rtmatch_min)))
      rtmatch_max <-
        .bincode(tdf[, "rtmax"], breaks = blim, include.lowest = TRUE)
      rtmatch_max <-
        ifelse(rtmatch_max < length(rttime),
               rtmatch_max + 1,
               rep(length(rttime), length(rtmatch_max)))
      rtmatch_med <-
        .bincode(tdf[, "rt"], breaks = blim, include.lowest = TRUE)

      # sel <- apply(rtmatch_min,rtmatch_max)

      tdf[, "rtmin"] <- rttime[rtmatch_min]
      tdf[, "rtmax"] <- rttime[rtmatch_max]
      tdf[, "rt"] <- rttime[rtmatch_med]

      ###only keep feature which have a different rt,rtmin and rtmax
      vdup <- apply(tdf[, c("rtmin", "rt", "rtmax")], 1, anyDuplicated)
      sel <- which(vdup == 0)

      tenv <- new.env()
      tenv[["chromPeaks"]] <- as.matrix(tdf[sel, ])
      mzdata@msFeatureData@.xData <- tenv
      mzdata@featureData@data <- mzraw@featureData@data
      mzdata@processingData@files <- path_raw
      return(list(mzdata, sel, sel_idx))
    }
    ldata <-
      convertToCliqueMS(dm,
                        idx = idx,
                        path_raw = raw_data,
                        mzraw = mzraw,
                        ref_xcms = ref_xcms)

    mzdata <- ldata[[1]]
    ###We just remove all the EICs value for simplicit


    ##We only extract EICs for feture which seems to fall in the correct boundaries
    sel <- ldata[[2]]
    sel_idx <- ldata[[3]]
    anclique <- createanClique(mzdata)


    createNetworkWithFilter <-
      function (mzdata,
                peaklist,
                filter = TRUE,
                mzerror = 5e-06,
                intdiff = 1e-04,
                rtdiff = 1e-04,
                cosFilter = 0.3)
      {
        sink("/dev/null")
        eicmat <- cliqueMS:::defineEIC(mzdata)
        sink(file = NULL)
        sparseeic <- as(t(eicmat), "sparseMatrix")
        cosTotal <- qlcMatrix::cosSparse(sparseeic)
        if (filter == TRUE) {
          filterOut <-
            cliqueMS:::filterFeatures(
              cosTotal,
              peaklist,
              mzerror = mzerror,
              rtdiff = rtdiff,
              intdiff = intdiff
            )
          cosTotal <- filterOut$cosTotal
          peaklist <- filterOut$peaklist
          # message(paste("Features filtered:", length(filterOut$deleted),
          #               sep = " "))
        }
        network <- igraph::graph.adjacency(cosTotal,
                                           weighted = TRUE,
                                           diag = FALSE,
                                           mode = "undirected")
        igraph::V(network)$id = seq_len(nrow(peaklist))
        nozeroEdges = igraph::E(network)[which(igraph::E(network)$weight >
                                                 cosFilter)]
        network <- igraph::subgraph.edges(network, nozeroEdges)
        igraph::E(network)$weight <- round(igraph::E(network)$weight,
                                           digits = 10)
        igraph::E(network)$weight[which(igraph::E(network)$weight ==
                                          1)] <- 0.99999999999
        return(list(network = network, peaklist = peaklist))
      }
    netlist <-
      createNetworkWithFilter(
        mzdata,
        anclique@peaklist[sel, ],
        filter = TRUE,
        mzerror = 1e-5,
        intdiff = 1e4,
        rtdiff = 1e-3,
        cosFilter = cosFilter
      )$network

    if (length(netlist) == 1 && is.na(netlist)) {
      return(matrix(nrow = 0, ncol = 3))
    }
    alle <- as_data_frame(netlist, "edges")
    vid <- vertex_attr(netlist, name = "id")
    alle[, 1] <- sel_idx[sel[vid[alle[, 1]]]]
    alle[, 2] <- sel_idx[sel[vid[alle[, 2]]]]
    alle <- as.matrix(alle)
    return(alle)
  }



createNetworkMultifiles <-
  function(dm,
           raw_files,
           opened_raw_files,
           match_files,
           size_batch = 5,
           ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData",
           cosFilter = 0.4,
           bpp = NULL) {
    size_batch <- min(size_batch, length(raw_files) - 1)
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
    seq_cut <- seq(1, length(raw_files), by = size_batch)
    if (seq_cut[length(seq_cut)] != length(raw_files)) {
      seq_cut <- c(seq_cut, length(raw_files) + 1)
    } else{
      seq_cut[length(seq_cut)] <- length(raw_files) + 1
    }
    # message("Building cosine similarity network", appendLF = FALSE)


    for (i in 1:(length(seq_cut) - 1)) {
      ###We compute the network for the selected files.

      ledges <-
          bpmapply(
            seq_cut[i]:(seq_cut[i + 1] - 1),
            as.list(raw_files[seq_cut[i]:(seq_cut[i + 1] - 1)]),
            as.list(opened_raw_files[seq_cut[i]:(seq_cut[i + 1] - 1)]),
            FUN = computeNetworkRawfile,
            MoreArgs = list(
              ref_xcms = ref_xcms,
              dm = dm,
              cosFilter = cosFilter
            ),
            BPPARAM = bpp
          #bptry( )
        )
      # message("ledges",format(object.size(ledges),"Mb"))


      found_errors <- sapply(ledges, function(x) {
        "remote_perror" %in% class(x)
      })
      ##Cahcing errors
      if (sum(found_errors)>=1){
        # browser()
        ledges <- attr(ledges, "result")

        # print(which(found_errors))
        ledges <-
          ledges[!found_errors]

        # if(length(found_errors)>1){
        #   print(sapply(ledges[found_errors],function(x,rr){attr(x,"traceback")},simplify = FALSE))
        # }
      }
      ###We add the edges to the main network
      for (j in seq_along(ledges)) {
        alle <- ledges[[j]]
        # message("Added ",nrow(alle),"edges")
        if (is.null(alle)) next
        if (nrow(alle) == 0)
          next
        alle <- as.matrix(alle)
        ##We update the cost matrix
        cosMat[alle[, c(1, 2)]] <-
          (cosMat[alle[, 1:2]] * countMat[alle[, 1:2]] + alle[, 3]) /
          (countMat[alle[, 1:2]] + 1)
        cosMat[alle[, c(2, 1)]] <- cosMat[alle[, c(1, 2)]]
        countMat[alle[, c(1, 2)]] <- countMat[alle[, c(1, 2)]] + 1
        # message("cosMat",format(object.size(cosMat),"Mb"))
      }
    }
    # message("\nDone")

    ##We get the first useless elemnts
    sel_val <- rowSums(cosMat) != 0
    # message(sum(sel_val),"signals with features.")
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
    # message("gadj",format(object.size(gadj),"Mb"))
    ldata <-
      convertToCliqueMS(dm,
                        path_raw = raw_files[[1]],
                        mzraw = opened_raw_files[[1]],
                        ref_xcms = ref_xcms)

    mzdata <- ldata[[1]]

    intensity <-
      apply(dm, 1, function(x) {
        mean(x[3:length(x)], na.rm = TRUE)
      })
    anclique <- createanClique(mzdata)
    anclique@peaklist <- data.frame(
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



    anclique@network <- gadj
    return(anclique)
  }


###Annotations of cliques
annotateCliques <- function(cliques,
                            peaklist,
                            adducts,
                            main_adducts,
                            ionization_mode,
                            val_int,
                            bpp = NULL) {
  if (is.null(bpp))
    bpp <- bpparam()
  for(ip in seq_along(cliques)){
    annotateCliqueInterpretMSspectrum(cliques[[ip]],
    dm = peaklist,
    adducts = adducts,
    main_adducts = main_adducts,
    ionization_mode = ionization_mode,
    val_int = val_int)
  }
  vfeat <-
    bplapply(
      cliques,
      FUN = annotateCliqueInterpretMSspectrum,
      dm = peaklist,
      adducts = adducts,
      main_adducts = main_adducts,
      ionization_mode = ionization_mode,
      val_int = val_int,
      BPPARAM = bpp
    )

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
  ###index is always the position of the feature in the originl data matrix.
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
  function(clique,
           adducts,
           main_adducts,
           ionization_mode,
           dm,
           val_int,
           ppm = 10,
           dmz = 0.005) {
    library(InterpretMSSpectrum)
    def_ion <- NULL
    if (ionization_mode == "positive") {
      def_ion <- "[M+H]+"
    } else{
      def_ion <- "[M-H]-"
    }
    if (length(clique) == 1) {
      return(list(
        data.frame(
          index = clique[1],
          mz = dm[clique[1], 1],
          intb = val_int[clique[1]],
          isogr = NA,
          iso = NA,
          charge = NA,
          adduct = def_ion,
          ppm = NA,
          label = def_ion
        )
      ))
    }
    sel_clust_idx <- clique
    sel_idx <- seq_along(sel_clust_idx)

    current_val <- cbind(dm[sel_clust_idx, 1], val_int[sel_clust_idx])
    colnames(current_val) <- c("mz", "int")

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
          adduct = def_ion,
          ppm = NA,
          label = def_ion
        )
        num_feat <- num_feat + 1
        break
      }
      
      annots <-
        tryCatch(suppressWarnings(findMAIN(
          current_val[sel_idx, , drop = FALSE],
          ionmode = ionization_mode,
          rules = adducts,
          adducthyp = main_adducts,
          mzabs = dmz,
          ppm = ppm,
          mainpkthr = 0.2
        ))[[1]],error=function(e){return(NA)})
      
      
      # val <- data.frame(mz=c(150,151.01,150-18,415),int=c(100,30,10,95))
      # findMAIN(val,ionmode="positive")
      
      # mz int isogr iso charge     adduct      ppm      label
      # 132.00  10    NA  NA     NA [M+H-H2O]+ 80.03561 [M+H-H2O]+
      # 150.00 100     1   0      1     [M+H]+       NA     [M+H]+
      # 151.01  30     1   1      1       <NA>       NA       <NA>
      # 415.00  95    NA  NA     NA       <NA>       NA       <NA>
      
      ####In case of error we consider all the values as the default
      
      
      
      sel_adducts <- which(!is.na(annots[, "adduct"]))
      if (length(sel_adducts) <= 0)
        break

      if (length(sel_adducts) == 1) {
        annots$label[sel_adducts] <- def_ion
        annots$adduct[sel_adducts] <- def_ion
      }
      
      c("index",
        "mz",
        "intb",
        "isogr",
        "iso",
        "charge",
        "adduct",
        "ppm",
        "label")
      
      
      

      sel_iso <- annots$isogr[sel_adducts]
      sel_iso <- sel_iso[!is.na(sel_iso)]

      sel_iso <- which(annots$isogr %in% sel_iso)
      sel_feat <- union(sel_adducts, sel_iso)

      all_features[[num_feat]] <-
        cbind(sel_clust_idx[sel_idx[sel_feat]], annots[sel_feat, ])
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

convertFeatures <- function(resAnnot, polarity = "negative") {
  ####We change all the data
  unique_groups <- unique(resAnnot$group_label)
  resList <- vector(mode = "list", length = length(unique_groups))
  current_dec <- 0
  def_label <- "[M+H]+"
  if (polarity == "negative") {
    def_label <- "[M-H]-"
  }
  message("Features conversion: ")
  for (ig in seq_along(resList)) {
    if (((ig * 10) %/% length(resList)) != current_dec) {
      message(current_dec * 10, " ", appendLF = FALSE)
      current_dec <- (ig * 10) %/% length(resList)
    }
    xannot <- resAnnot[resAnnot$group_label == unique_groups[ig], ]
    vannot <- unique(xannot[, "isogr"])
    vannot <- vannot[!is.na(vannot)]
    p_unique <- NULL
    if (length(vannot) != 0) {
      p_unique <- match(vannot, xannot[, "isogr"])
    } else{
      p_unique <- match(vannot, xannot[, "isogr"])
    }

    sel_adducts <- xannot[, "label"]
    vvv <- which(!is.na(xannot[, "adduct"]))
    pmax <- vvv[which.max(xannot[vvv, "intb"])]
    vs <-
      InterpretMSSpectrum:::getRuleFromIonSymbol(xannot[pmax, "label"])
    neutral_mass <-
      as.numeric((vs[2] * xannot[pmax, "mz"] - vs[4]) / abs(vs[3]))
    tempdf <-
      data.frame(
        index = xannot[, "index"],
        neutral_mass = rep(neutral_mass, nrow(xannot)),
        label = as.character(xannot[, "label"]),
        charge = xannot[, "charge"],
        group_label = xannot[, "group_label"],
        adduct = as.character(xannot[, "label"]),
        ref_feature = rep(xannot[pmax, "index"], nrow(xannot)),
        stringsAsFactors = FALSE
      )
    if (length(vannot) == 0) {
      resList[[ig]] <- tempdf
      next
    }
    for (p in seq_along(p_unique)) {
      pp <- p_unique[p]
      pisos <- which(xannot[, "isogr"] == xannot[pp, "isogr"])
      ref_label <- as.character(xannot[pp, "label"])
      if (length(pisos) == 1) {
        tempdf[pisos, "adduct"] <- as.character(ref_label)
        next
      }
      if (is.na(ref_label))
        ref_label <- def_label

      ###We update the label of all the peaks
      labelC13 <-
        paste(ref_label, c("", paste("+", 1:(
          length(pisos) - 1
        ), sep = "")))
      tempdf[pisos, "adduct"] <- labelC13
      tempdf[pisos, "label"] <- labelC13
    }
    ###We reorder the data.frame
    tempdf <- tempdf[order(tempdf$index), ]
    resList[[ig]] <- tempdf
  }
  return(resList)
}


###The full data matrix annotated.
buildDataMatrixFull <- function(dm, annot, path_output,
                                num_line = 15000) {
  ###We write by batch to avoid memory oerloads.
  cnames <-
    c("mz",
      "rt",
      "main_peak",
      "annotations",
      "neutral_mass",
      cnames[3:ncol(dm)])
  f <- file(path_output, "a")

  seq_line <- seq(1, length(annot), by = num_line)
  if (seq_line[length(seq_line)] != length(annot)) {
    seq_line <- c(seq_line, length(annot))
  } else{
    seq_line[length(seq_line)] <- length(annot) + 1
  }
  message("Building full data-matrix in ",
          length(seq_line) - 1,
          " batches.")
  write_col <- TRUE

  for (i in 1:(length(seq_line) - 1)) {
    sel_idx <- seq_line[i]:(seq_line[i + 1] - 1)
    vannot <- sapply(annot[sel_idx], function(an, dm, cnames) {
      sub_dm <-
        cbind(dm[an[, "index"], c(1, 2)], an[, "adduct"], an[, "group_label"],
              an[, "neutral_mass"], dm[an[, "index"], 3:ncol(dm)])
      colnames(sub_dm) <-
        c("mz",
          "rt",
          "annotation",
          "group",
          "neutral_mass",
          cnames[3:ncol(dm)])
      return(sub_dm)
    }, dm = dm, cnames = colnames(dm),
    simplify = FALSE)
    vannot <- do.call(rbind, vannot)
    write.table(
      vannot,
      f,
      sep = ";",
      col.names = write_col,
      row.names = FALSE
    )
    write_col <- FALSE
  }
  ###Write the full dtaa matrix to files sequentially.
  invisible(close(f))
}


buildDataMatrixSimplified <-
  function(dm,
           annot,
           path_output,
           vsep = ";",
           num_line = 15000) {
    ###We write by batch to avoid memory oerloads.
    cnames <-
      c("mz",
        "rt",
        "main_peak",
        "annotations",
        "neutral_mass",
        cnames[3:ncol(dm)])
    f <- file(path_output, "a")


    seq_line <- seq(1, length(annot), by = num_line)
    if (seq_line[length(seq_line)] != length(annot)) {
      seq_line <- c(seq_line, length(annot))
    } else{
      seq_line[length(seq_line)] <- length(annot) + 1
    }
    message("Building simplified data-matrix in ",
            length(seq_line) - 1,
            " batches.")
    write_col <- TRUE
    ### This is done for every batch
    for (i in 1:(length(seq_line) - 1)) {
      sel_idx <- seq_line[i]:(seq_line[i + 1] - 1)
      vannot <- sapply(annot[sel_idx], function(an, dm, cnames) {
        refv <- an[1, "ref_feature"]
        main_adduct <- an[match(refv, an[, 1]), "adduct"]
        all_adducts <- paste(an[, "adduct"], collapse = "|")
        neutral_mass <- an[1, "neutral_mass"]
        tdf <-
          data.frame(
            mz = dm[refv, 1],
            rt = dm[refv, 2],
            main_adduct = main_adduct,
            all_adducts = all_adducts,
            neutral_mass = neutral_mass
          )
        tdf <- cbind(tdf, dm[refv, 3:ncol(dm)])
        colnames(tdf) <-
          c("mz",
            "rt",
            "main_peak",
            "annotations",
            "neutral_mass",
            cnames[3:ncol(dm)])
        return(tdf)
      }, dm = dm, cnames = colnames(dm),
      simplify = FALSE)
      vannot <- do.call(rbind, vannot)
      write.table(
        vannot,
        f,
        sep = ";",
        col.names = write_col,
        row.names = FALSE
      )
      write_col <- FALSE
    }
    invisible(close(f))
  }


groupFeatures <-
  function(dm,
           val_int,
           raw_files,
           opened_raw_files,
           match_files,
           adducts,
           main_adducts,
           ionization_mode,
           ppm = 10,
           dmz = 0.005,
           size_batch = 10,
           cut_size = 10000,
           cosFilter = 0.5,
           ref_xcms = "X:/Documents/dev/script/diverse/xcms_raw_with_peaks.RData",
           path_matching = "cliques_matching.cpp",
           bbp = NULL) {
    if (is.null(bpp))
      bpp <- bpparam()

    #We have to source the matchign function.
    sourceCpp(path_matching)

    ###We sort dm by retnetion time
    ort <- order(dm[, "rt"], decreasing = FALSE)


    ###We plit the files evnetually.
    number_features <- nrow(dm)
    cut_rts <- seq(1, number_features, by = cut_size/2)
    if (cut_rts[length(cut_rts)] != number_features) {
      cut_rts <- c(cut_rts, number_features + 1)
    } else{
      cut_rts[length(cut_rts)] <- cut_rts[length(cut_rts)] + 1
    }
    if(length(cut_rts)==2) cut_rts <- c(cut_rts,cut_rts[length(cut_rts)])

    ###WE update the cliaues vector at every step
    cliques <- vector(mode = "list", length = 10000)
    assignments <- rep(NA_integer_, nrow(dm))
    size <- rep(0L, nrow(dm))
    ###The first id needs to be -1
    current_id <- as.integer(-1)
    ##We update all the cliques.
    message("Annotations task divided in ",length(cut_rts)-2," batches.")
    for (i in seq(1, length(cut_rts) - 2)) {
      message("Processing batch ",i)
      # message("Variables:",cut_rts[i],"-",cut_rts[i+2],"current cliques size is ",length(cliques))
      sel_idx <- ort[cut_rts[i]:cut_rts[i + 2]]
      anclique <- createNetworkMultifiles(
        dm[sel_idx, , drop = FALSE],
        raw_files,
        opened_raw_files,
        match_files,
        size_batch = size_batch,
        ref_xcms = ref_xcms,
        cosFilter = cosFilter,
        bpp = bpp
      )
      ###we compute the cliques
      # message("Computing cliques")
      sink(file="/dev/null")
      anclique <- computeCliques(anclique, 1e-5, TRUE)
      sink(NULL)
      ###We correct the index for subselection.
      for (ic in seq_along(anclique@cliques)) {
        anclique@cliques[[ic]] <- sel_idx[anclique@cliques[[ic]]]
      }
      # message("Merging cliques")
      
      ###TO DEBUG only
      # res_list <- list(cliques, anclique@cliques, assignments, size, current_id)
      # out_path <- file.path(Sys.getenv('OUTPUT'),"save_merge.rds")
      # saveRDS(res_list,file = out_path)
      
      # ocliques <- cliques
      cliques <-
        mergeCliques(cliques, anclique@cliques, assignments, size, current_id)

    }

    ###We convert the cliques ID back to the unsorted point
    # for (ic in seq_along(cliques)) {
    #   cliques[[ic]] <- cliques[[ic]]
    # }

    summarized_df <- data.frame(
      mz = dm[, "mz"],
      mzmin = dm[, "mz"] - (dm[, "max_mz"] + dm[, "min_mz"]) * 1.5 ,
      mzmax = dm[, "mz"] + (dm[, "max_mz"] - dm[, "min_mz"]) * 1.5,
      rt = 60 * (dm[, "rt"]),
      rtmin = 60 * (dm[, "rt"] - dm[, "mean_peakwidth"] - 0.002),
      rtmax = 60 * (dm[, "rt"] + dm[, "mean_peakwidth"] + 0.002),
      into = val_int,
      intb = val_int,
      maxo = val_int,
      sn = rep(10, nrow(dm)),
      sample =
        rep(1, nrow(dm)),
      is_filled = rep(0, nrow(dm))
    )
    ###We just have to construct the definitive peaktable.
    ###We can redirect all the cliques in the data
    empty_graph <- function() {
      temp <- matrix(rnorm(4), nrow = 2)
      ga <- graph_from_adjacency_matrix(adjmatrix = temp,
                                        mode = "undirected",
                                        weighted = TRUE)
      return(ga)
    }

    anclique@network <- empty_graph()
    pint <- getIntensityPos(dm)[1]

    message("Annotating")
    res_df <- annotateCliques(cliques,
                              summarized_df,
                              adducts,
                              main_adducts,
                              ionization_mode,
                              val_int,
                              bpp = bpp)

    ####We add the annotation for all the peaks
    annot <- convertFeatures(res_df, polarity = ionization_mode)
    message("Converting features.")
    return(annot)
  }

####Actual processing in the pipeline.
args <- commandArgs(trailingOnly = TRUE)
# args <- c("U:/users/Alexis/data/all_2mins/res_neg_2/datamatrices/datamatrix_84b6673cb2b08293f0365e15464a88ed.csv",
#           "U:/users/Alexis/data/all_2mins/res_neg_2/processing_db.sqlite",
#           "/U:/users/Alexis/data/all_2mins/res_neg_2/datamatrices/annotated_peaktable_84b6673cb2b08293f0365e15464a88ed_full.csv",
#           "U:/users/Alexis/data/all_2mins/res_neg_2/datamatrices/annotated_peaktable_84b6673cb2b08293f0365e15464a88ed_reduced.csv",
#           "5",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/xcms_raw_model.RData",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/adducts_neg.txt",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/adducts_main_neg.txt",
#           "negative",
#           "15.0",
#           "0.007",
#           "50",
#           "2",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/Rscript/cliques_matching.cpp")



# print(args)
# args <- c("U:/users/Alexis/examples_lcms_workflow/output/datamatrices/datamatrix_3063282bcc8e77018d0b6912579a4115.csv",
#           "U:/users/Alexis/examples_lcms_workflow/output/processing_db.sqlite",
#           "E:/out_full_d.csv",
#           "E:/out_simple_d.csv",
#           "5",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/xcms_raw_model.RData",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/adducts_pos.txt",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/adducts_main_pos.txt",
#           "positive",
#           "15",
#           "0.01",
#           "50",
#           "2",
#           "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/Rscript/cliques_matching.cpp")

PATH_DATAMATRIX <- args[1]
PATH_DB <- args[2]
PATH_OUTPUT_FULL <- args[3]
PATH_OUTPUT_SIMPLE <- args[4]
NUM_CORES <- ceiling(as.numeric(args[5]) / 2)
PATH_MODEL <- args[6]
PATH_ADDUCTS <- args[7]
PATH_MAIN_ADDUCTS <- args[8]
POLARITY <- args[9]
PPM <-  as.numeric(args[10])
DMZ <-  as.numeric(args[11])
###We peak the FILE_USED most intense files.
FILES_USED <- min(as.numeric(args[12]), 25)
FILTER_NUMS <- max(1, as.numeric(args[13]))
PATH_MATCHING <- args[14]
NUM_BY_BATCH <- 7000
if (length(args) == 15) {
  NUM_BY_BATCH <- as.numeric(args[15])
}
#

###Debuggiong the test on karin data
# PATH_DATAMATRIX <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res_thin/datamatrices/datamatrix_b01904e0a545f19857549535352af3b9.csv"
# PATH_DB <- "U:/users/Alexis/data/data_karin/mml_v2/neg/res/processing_db.sqlite"
# PATH_OUTPUT_FULL <- "E:/out_full.csv"
# PATH_OUTPUT_SIMPLE <- "E:/out_simple.csv"
# NUM_CORES <- as.numeric(4)
# PATH_MODEL <- "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/xcms_raw_model.RData"
# PATH_ADDUCTS <- "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/adducts_neg.txt"
# PATH_MAIN_ADDUCTS <- "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/data/adducts_main_neg.txt"
# POLARITY <- "negative"
# PPM <-  as.numeric(20)
# DMZ <-  as.numeric(0.02)
# FILTER_NUMS <- max(1,as.numeric(2))
# FILES_USED <- 4
# PATH_MATCHING <- "C:/Users/dalexis/Documents/python/lcmsprocessing/pylcmsprocessing/Rscript/cliques_matching.cpp"
# NUM_BY_BATCH <- 20000


##Debugging standard workflow

##reading data matrices
dm <- read.table(PATH_DATAMATRIX, header = TRUE, sep = ",")

posIntensities <- getIntensityPos(dm)
num_detect <-
  apply(dm[, posIntensities, drop = FALSE], 1, function(x) {
    sum(!is.na(x))
  })

vdetect <- num_detect >= FILTER_NUMS
###We only keep the fitting ammount of features.
while (sum(vdetect) > 200000) {
  FILTER_NUMS <- FILTER_NUMS + 1
  vdetect <- num_detect >= FILTER_NUMS
}

# message("Retained ", sum(vdetect), " signals on ", nrow(dm))
dm <- dm[vdetect, , drop = FALSE]

###Reading the raw files
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
raw_files <- dbGetQuery(dbb, "SELECT path FROM samples")[, 1]
dbDisconnect(dbb)
# raw_files <- str_replace(raw_files,pattern = "/sauer1",replacement = "U:")


####Selecting the msot intense files
val_int <- apply(dm[, posIntensities], 2, sum, na.rm = TRUE)
sel_files <-
  order(val_int, decreasing = TRUE)[1:min(FILES_USED, length(val_int))]
raw_files <- raw_files[sel_files]


val_int_var <- apply(dm[, posIntensities], 1, mean, na.rm = TRUE)
###Setting up the parallel processing
bpp <- NULL
if (get_os() == "win") {
  bpp <- SnowParam(workers = NUM_CORES)
} else{
  bpp <- MulticoreParam(workers = min(NUM_CORES, 4))
}

opened_raw_files <- sapply(raw_files,readMSData, mode = "onDisk")

fadd <- file(PATH_ADDUCTS, "r")
adducts <- readLines(fadd)
close(fadd)

###
fadd <- file(PATH_MAIN_ADDUCTS, "r")
main_adducts <- readLines(fadd)
main_adducts <- main_adducts[sapply(main_adducts,startsWith,prefix="[M")]
###The double molecule are removed
close(fadd)

####We map the sample name of the vector data
base_sample <-
  str_split(basename(raw_files),
            pattern = fixed("\\."),
            simplify = TRUE)[, 1]

##Give the posaition of raw_file on the data matrix
cnames <- colnames(dm)
cnames <- str_sub(cnames, 5, -5)
vav <- adist(base_sample, cnames)
match_files <- apply(vav, 1, which.min)

# match_files <- sapply(base_sample,function(x,vref){which(endsWith(vref,suffix=x))},vref=colnames(dm),simplify=TRUE)


###If the software is crashing we divied the number of feature by 2 eventually
# print(head(raw_files))
annot <-
  groupFeatures(
    dm[, c("mz", "rt","min_mz","max_mz","min_rt","max_rt","mean_peakwidth",colnames(dm)[posIntensities[sel_files]])],
    val_int_var,
    raw_files,
    opened_raw_files,
    match_files,
    adducts,
    main_adducts,
    ionization_mode = POLARITY,
    ppm = PPM,
    dmz = DMZ,
    size_batch = NUM_CORES,
    cut_size = NUM_BY_BATCH,
    cosFilter = 0.6,
    ref_xcms = PATH_MODEL,
    path_matching = PATH_MATCHING,
    bbp = bpp
  )

###If the file already exists at this step we erase it
if(file.exists(PATH_OUTPUT_SIMPLE)) file.remove(PATH_OUTPUT_SIMPLE)
dm_full <- buildDataMatrixSimplified(dm, annot, PATH_OUTPUT_SIMPLE)
if(file.exists(PATH_OUTPUT_FULL)) file.remove(PATH_OUTPUT_FULL)
dm_full <- buildDataMatrixFull(dm, annot, PATH_OUTPUT_FULL)
