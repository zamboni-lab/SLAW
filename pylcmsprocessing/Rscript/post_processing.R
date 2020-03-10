###Given as input a list of file or a list of ids, extract the associated EICs.
suppressWarnings(suppressMessages(library(DBI, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(RSQLite, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(stringr, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(igraph, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))

args <- commandArgs(trailingOnly = TRUE)
# print(args)
# args <- c(
#   "/sauer1/users/Alexis/examples_lcms_workflow/input/target.csv",
#   "/sauer1/users/Alexis/examples_lcms_workflow/output/temp/raw_files.txt",
#   "/sauer1/users/Alexis/examples_lcms_workflow/output/processing_db.sqlite",
#   "/sauer1/users/Alexis/examples_lcms_workflow/output/hdf5/eics.hdf53063282bcc8e77018d0b6912579a4115.pdf",
#   "/sauer1/users/Alexis/examples_lcms_workflow/output/figures/peaks3063282bcc8e77018d0b6912579a4115.pdf",
#   "/sauer1/users/Alexis/examples_lcms_workflow/output/figures/diagnosis3063282bcc8e77018d0b6912579a4115.pdf",
#   "/sauer1/users/Alexis/examples_lcms_workflow/output/datamatrices/targetted_rt3063282bcc8e77018d0b6912579a4115.csv",
#   "/sauer1/users/Alexis/examples_lcms_workflow/output/datamatrices/targetted_int3063282bcc8e77018d0b6912579a4115.csv",
#   "0.05",
#   "0.03",
#   "5"
# )

# args <- str_replace(args,"/sauer1","U:")




INPUT_FEATURES <- args[1]
##For later use
INPUT_RAW_FILES <- args[2]
PATH_DB <- args[3]
###For later use.
OUTPUT_HDF5 <- args[4]
OUTPUT_TARGET_PDF <- args[5]
OUTPUT_SUMMARY_PDF <- args[6]
OUTPUT_TARGETTED_RT_TABLE <- args[7]
OUTPUT_TARGETTED_INT_TABLE <- args[8]
TOL_MZ <- as.numeric(args[9])
TOL_RT <- as.numeric(args[10])
NCORES <- as.numeric(args[11])

MARGIN_RT <- 2 / 60 #(3s)
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

###We read the assocaited data.frame
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
PATH_DATAMATRIX <-
  dbGetQuery(dbb, "SELECT annotated_peaktable_reduced FROM PEAKPICKING")[1, 1]
dbDisconnect(dbb)


dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
raw_files <- dbGetQuery(dbb, "SELECT path FROM samples WHERE level='MS1'")[, 1]
dbDisconnect(dbb)

dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
peaktables <-
  dbGetQuery(dbb, "SELECT output_ms FROM processing WHERE output_ms!='NOT PROCESSED'")[, 1]
dbDisconnect(dbb)
#
# peaktables <- str_replace(peaktables,"/sauer1","U:")
# raw_files <- str_replace(raw_files,"/sauer1","U:")
# PATH_DATAMATRIX <- str_replace(PATH_DATAMATRIX,"/sauer1","U:")

if (length(raw_files) > 100) {
  sint <- sample.int(length(raw_files), size = 100)
  raw_files <- raw_files[sint]
  peaktables <- peaktables[sint]
}

traw_files <- str_split(raw_files, fixed("."), simplify = TRUE)
traw_files <- basename(traw_files[, 1])

bpp <- NULL
if (get_os() == "win") {
  bpp <- SnowParam(workers = NCORES)
} else{
  bpp <- MulticoreParam(workers = min(NCORES, 4))
}

##We know extract the value
summarize <- function(pt, val) {
  library(data.table)
  pt <- fread(pt, header = TRUE, sep = ",")
  apply(pt[, ..val, drop = FALSE], 2, function(x) {
    sel <- !is.na(x) & !is.infinite(x) & !x > 1e11
    mean(x[sel])
  })
}


sum_metrics <-
  bplapply(
    split(peaktables, f = seq_along(peaktables)),
    summarize,
    val = c("SN", "intensity", "peakwidth", "right_on_left_assymetry"),
    BPPARAM = bpp
  )


suppressWarnings(suppressMessages(library(ropls, warn.conflicts = FALSE)))
sum_metrics <- do.call(rbind, sum_metrics)

sel_unique <- !duplicated(sum_metrics)

###Stupid case of 2 times the same sample.
sum_metrics <- sum_metrics[sel_unique,]

row.names(sum_metrics) <- traw_files[sel_unique]
sum_metrics_scaled <- scale(sum_metrics)
pdf(OUTPUT_SUMMARY_PDF)
sink("/dev/null")
val <- opls(
  sum_metrics_scaled,
  y = NULL,
  crossvalI = 1,
  predI = 2
)
plot(
  sum_metrics[, c("peakwidth", "right_on_left_assymetry")],
  xlab = "Mean peakwidth",
  ylab = "Mean assymetry",
  main = "Peakshape diagnosis",
  pch = 16,
  cex = 0.8,
  col = "darkred"
)
text(
  x = sum_metrics[, "peakwidth"],
  y = sum_metrics[, "right_on_left_assymetry"],
  labels = row.names(sum_metrics),
  col = "darkred",
  cex = 0.7,
  adj = c(0, 0)
)
plot(
  sum_metrics[, c("SN", "intensity")],
  xlab = "Mean SN",
  ylab = "Mean intensity",
  main = "Intensity diagnosis",
  pch = 16,
  cex = 0.8,
  col = "darkblue"
)
text(
  x = sum_metrics[, "SN"],
  y = sum_metrics[, "intensity"],
  labels = row.names(sum_metrics),
  col = "darkblue",
  cex = 0.7,
  adj = c(0, 0)
)
dev.off()
sink(file = NULL)




matchLCMSsignals <-
  function(mz_data,
           rt_data,
           mz_ref,
           rt_ref,
           tol_mz = 0.01,
           tol_rt = 10) {
    ###We map the function to the napping
    # vm <- mineMS2:::matchMzs(mz_ref, mz_data, ppm = ppm,dmz=dmz)
    vm <- xcms:::fastMatch(mz_ref, mz_data, tol = tol_mz)
    pf <-  which(!sapply(vm, is.null))

    createEdge <- function(x,
                           xi,
                           yi,
                           mz_data,
                           rt_data,
                           mz_ref,
                           rt_ref,
                           alpha_rt,
                           intlog) {
      if (is.na(x)) {
        return(
          data.frame(
            from = numeric(0),
            to = numeric(0),
            weight = numeric(0),
            mz1 = numeric(0),
            rt1 = numeric(0),
            mz2 = numeric(0),
            rt2 = numeric(0)
          )
        )
      } else{
        ###deviation in ppm
        cmz <- abs(mz_ref[xi] - mz_data[x]) / tol_mz
        crt <- abs(rt_ref[xi] - rt_data[x])
        if (all(crt > (2 * tol_rt)))
          return(
            data.frame(
              from = numeric(0),
              to = numeric(0),
              weight = numeric(0),
              mz1 = numeric(0),
              rt1 = numeric(0),
              mz2 = numeric(0),
              rt2 = numeric(0)
            )
          )

        cmz <- max(cmz) * 1.5 - cmz
        crt <- max(crt) * 1.1 - crt

        # cint <-  log10(int_data[x]) / intlog
        score <- cmz + crt# + cint


        return(
          data.frame(
            from = rep(xi, length(x)),
            to = x + length(mz_ref),
            weight = score,
            mz1 = mz_ref[rep(xi, length(x))],
            rt1 = rt_ref[rep(xi, length(x))],
            mz2 = mz_data[x],
            rt2 = rt_data[x]
          )
        )
      }
    }

    df_edges <-
      mapply(
        vm[pf],
        pf,
        seq_along(pf),
        FUN = createEdge,
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
      make_bipartite_graph(
        types = c(rep(FALSE, length(mz_ref)), rep(TRUE, length(mz_data))),
        edges = as.numeric(t(as.matrix(df_edges[, c(1, 2)]))),
        directed = FALSE
      )

    E(ge)$weight <- df_edges[, 3]

    ###IT WORKS
    mm <- igraph::max_bipartite_match(ge)
    ###We just add the matching part
    res_matching <-  mm$matching[1:length(mz_ref)]

    res_matching <-  sapply(res_matching, function(x, shift) {
      if (is.na(x))
        return(x)
      return(x - shift)
    }, shift = length(mz_ref))
    ###Output : a vector giving the postion of mz_ref in the data
    return(res_matching)
  }


matchMSsignals <-
  function(mz_data,
           mz_ref,
           num_detect,
           tol_mz = 0.01,
           tol_rt = 10) {
    ###We map the function to the napping
    # vm <- mineMS2:::matchMzs(mz_ref, mz_data, ppm = ppm,dmz=dmz)
    vm <- xcms:::fastMatch(mz_ref, mz_data, tol = tol_mz)
    found <- !sapply(vm, is.null)


    vfound <-
      mapply(
        vm[found],
        mz_ref[found],
        FUN = function(x, y, vmz, dd) {
          x[which.max(log10(dd[x]) * abs(vmz[x] - y))]
        },
        MoreArgs = list(vmz = mz_data, dd = num_detect)
      )


    res <- rep(NA, length(mz_ref))
    res[found] <- vfound
    return(res)

  }





vmz <- NULL
vrt <- NULL

dm <- read.table(PATH_DATAMATRIX, sep = ",", header = TRUE)
matched_features <- NULL
found_signals <- NULL
vfound <- NULL
if (file.exists(INPUT_FEATURES)) {
  input_signals <-
    read.table(
      INPUT_FEATURES,
      sep = ";",
      header = TRUE,
      stringsAsFactors = FALSE
    )
  if (ncol(input_signals) == 1) {
    ###Then we consider that the data have been directly removed.
    sel_idx <- input_signals[, 1]
    # vpp <- log10(apply(dm[,22:ncol(dm)],1,mean,na.rm=TRUE))
    # sel_idx <- sample.int(nrow(dm),10,prob = vpp/sum(vpp))
    vmz <- dm[sel_idx, "mz"]
    vrt <- dm[sel_idx, "rt"]
  } else if (ncol(input_signals) >= 2) {
    vmz <- input_signals[, "neutral_mass"]
    matched_features <- NULL
    if (!"rt" %in% colnames(input_signals)) {
      matched_features <- matchMSsignals(dm[, "neutral_mass"],
                                         vmz, TOL_MZ,
                                         num_detect = dm[, "num_detection"])
    } else{
      matched_features <- matchLCMSsignals(dm[, "neutral_mass"],
                                           dm[, "rt"],
                                           input_signals[, "neutral_mass"], input_signals[, "rt"], TOL_MZ,
                                           TOL_RT)
    }
    ###In very case we map the feature to the data.
    vmz <- dm[matched_features, "mz"]
    vrt <- dm[matched_features, "rt"]
  }

  found_signals <- !is.na(vmz)
  if(sum(found_signals)==0) stop("No metabolites found, ")

  titles <- rep("", length(vmz))
  if (ncol(input_signals) >= 2) {
    if ("name" %in% colnames(input_signals)) {
      titles <- input_signals[, "name"]
    } else{
      titles <- input_signals[, 3]
    }
  }

  supp_titles <-
    paste(
      "RT: ",
      sprintf("%0.2f", vrt - MARGIN_RT),
      "-",
      sprintf("%0.2f", vrt + MARGIN_RT),
      "\nMZ:",
      sprintf("%0.4f", vmz - MARGIN_MZ),
      "-",
      sprintf("%0.4f", vmz + MARGIN_MZ),
      sep = ""
    )
  titles <- paste(titles, supp_titles, sep = "\n")




  ####
  extract_EIC_raw_file <-
    function(path_xraw,
             path_peaktable,
             mz,
             matched = NULL,
             rt = NULL,
             mz_margin,
             rt_margin,
             num_detect = NULL) {
      ###We also map the peaktbale in the data
      suppressWarnings(suppressMessages(library(
        xcms, quietly = TRUE, warn.conflicts = FALSE
      )))
      suppressWarnings(suppressMessages(library(
        igraph, quietly = TRUE, warn.conflicts = FALSE
      )))

      if (!file.exists(path_peaktable))
        stop(paste(path_peaktable, "file does not exist"))
      peaktable <-
        read.table(path_peaktable, header = TRUE, sep = ",")
      if (!file.exists(path_xraw))
        stop(paste(path_xraw, "file does not exist"))
      xraw <- suppressMessages(suppressWarnings(xcmsRaw(path_xraw)))

      matchLCMSsignals <-
        function(mz_data,
                 rt_data,
                 mz_ref,
                 rt_ref,
                 tol_mz = 0.01,
                 tol_rt = 10) {
          ###We map the function to the napping
          # vm <- mineMS2:::matchMzs(mz_ref, mz_data, ppm = ppm,dmz=dmz)
          vm <- xcms:::fastMatch(mz_ref, mz_data, tol = tol_mz)
          pf <-  which(!sapply(vm, is.null))

          createEdge <- function(x,
                                 xi,
                                 yi,
                                 mz_data,
                                 rt_data,
                                 mz_ref,
                                 rt_ref,
                                 alpha_rt,
                                 intlog) {
            if (is.na(x)) {
              return(
                data.frame(
                  from = numeric(0),
                  to = numeric(0),
                  weight = numeric(0),
                  mz1 = numeric(0),
                  rt1 = numeric(0),
                  mz2 = numeric(0),
                  rt2 = numeric(0)
                )
              )
            } else{
              ###deviation in ppm
              cmz <- abs(mz_ref[xi] - mz_data[x]) / tol_mz
              crt <- abs(rt_ref[xi] - rt_data[x])
              if (all(crt > (2 * tol_rt)))
                return(
                  data.frame(
                    from = numeric(0),
                    to = numeric(0),
                    weight = numeric(0),
                    mz1 = numeric(0),
                    rt1 = numeric(0),
                    mz2 = numeric(0),
                    rt2 = numeric(0)
                  )
                )

              cmz <- max(cmz) * 1.5 - cmz
              crt <- max(crt) * 1.1 - crt

              # cint <-  log10(int_data[x]) / intlog
              score <- cmz + crt# + cint


              return(
                data.frame(
                  from = rep(xi, length(x)),
                  to = x + length(mz_ref),
                  weight = score,
                  mz1 = mz_ref[rep(xi, length(x))],
                  rt1 = rt_ref[rep(xi, length(x))],
                  mz2 = mz_data[x],
                  rt2 = rt_data[x]
                )
              )
            }

            # "from","to","weight"

          }

          df_edges <-
            mapply(
              vm[pf],
              pf,
              seq_along(pf),
              FUN = createEdge,
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
            make_bipartite_graph(
              types = c(rep(FALSE, length(mz_ref)), rep(TRUE, length(mz_data))),
              edges = as.numeric(t(as.matrix(df_edges[, c(1, 2)]))),
              directed = FALSE
            )

          E(ge)$weight <- df_edges[, 3]

          ###IT WORKS
          mm <- igraph::max_bipartite_match(ge)
          ###We just add the matching part
          res_matching <-  mm$matching[1:length(mz_ref)]

          res_matching <-  sapply(res_matching, function(x, shift) {
            if (is.na(x))
              return(x)
            return(x - shift)
          }, shift = length(mz_ref))
          ###Output : a vector giving the postion of mz_ref in the data
          return(res_matching)
        }


      matchMSsignals <-
        function(mz_data,
                 mz_ref,
                 num_detect,
                 tol_mz = 0.01,
                 tol_rt = 10) {
          ###We map the function to the napping
          # vm <- mineMS2:::matchMzs(mz_ref, mz_data, ppm = ppm,dmz=dmz)
          vm <- xcms:::fastMatch(mz_ref, mz_data, tol = tol_mz)
          found <- !sapply(vm, is.null)


          vfound <-
            mapply(
              vm[found],
              mz_ref[found],
              FUN = function(x, y, vmz, dd) {
                x[which.max(log10(dd[x]) * abs(vmz[x] - y))]
              },
              MoreArgs = list(vmz = mz_data, dd = num_detect)
            )


          res <- rep(NA, length(mz_ref))
          res[found] <- vfound
          return(res)

        }

      matched_features <- NULL
      if (missing(rt)) {
        matched_features <- matchMSsignals(peaktable[, 1],
                                           mz, mz_margin, num_detect =
                                             num_detect)
      } else{
        matched_features <- matchLCMSsignals(peaktable[, 1],
                                             peaktable[, 2],
                                             mz, rt, mz_margin,
                                             rt_margin)
      }


      matched <- !is.na(matched_features)

      #
      # ccc <- cbind(peaktable[matched_features[matched],1],
      #       peaktable[matched_features[matched],2],
      #       mz[matched],rt[matched])
      # print(apply(ccc,1,paste,collapse="_"))

      ###We select the feature which are matched
      ###we build a feature list.

      rlmzr <-
        split(cbind(peaktable[matched_features[matched], "mz_min"],
                    peaktable[matched_features[matched], "mz_max"]), f =
                1:sum(matched))
      rlrtr <-
        split(cbind((peaktable[matched_features[matched], "rt_min"]) * 60,
                    (peaktable[matched_features[matched], "rt_max"]) *
                      60), f = 1:sum(matched))

      lmzr <-
        split(cbind(peaktable[matched_features[matched], "mz_min"] - mz_margin,
                    peaktable[matched_features[matched], "mz_max"] + mz_margin),
              f = 1:sum(matched))
      lrtr <-
        split(cbind((peaktable[matched_features[matched], "rt_min"] - rt_margin) *
                      60,
                    (peaktable[matched_features[matched], "rt_max"] +
                       rt_margin) * 60), f = 1:sum(matched))
      xraw <- suppressMessages(suppressWarnings(xcmsRaw(path_xraw)))
      reslist <- vector(mode = "list", length = length(mz))

      vres <-
        mapply(
          lmzr,
          lrtr,
          rlmzr,
          rlrtr,
          FUN = function(mzr, rtr, rmzr, rrtr, xraw) {
            reic <- rawEIC(xraw, mzrange = mzr, rtrange = rtr)
            ###We compute the limit in scan.
            integ <- ifelse(((xraw@scantime[reic$scan] <= rrtr[2]) &
                               (xraw@scantime[reic$scan] >= rrtr[1])), rep(TRUE, length(reic$scan)),
                            rep(FALSE, length(reic$scan)))

            return(list(
              time = xraw@scantime[reic$scan],
              intensity = reic$intensity,
              peak = integ
            ))
          },
          MoreArgs = list(xraw = xraw),
          SIMPLIFY = FALSE
        )

      reslist[matched] <- vres
      ### WE return the limits of integration
      raw_rt <- rep(NA,length(matched_features))
      raw_rt[matched] <- peaktable[matched_features[matched], 2]

      return(list(reslist,raw_rt))
    }

  ###We now change the vilsualization of the software by team.



  plotPeaks <-
    function(eics, titles, names_raw) {
      #,name_files,name_compounds

      colors <- rainbow(length(eics))

      for (i in seq_along(eics[[1]])) {
        rtlim <- sapply(eics, function(x, idx) {
          if (is.null(x[[idx]]))
            return(c(NA_real_, NA_real_, NA_real_))
          c(range(x[[idx]]$time), max(x[[idx]]$intensity))
        }, idx = i)
        rtmin <- suppressWarnings(min(rtlim[1,], na.rm = TRUE))
        rtmax <- suppressWarnings(max(rtlim[2,], na.rm = TRUE))
        intmax <- suppressWarnings(max(rtlim[3,], na.rm = TRUE))
        if (is.infinite(rtmin))
          next

        ####We build the legend vector
        sel_raw <- which(!is.na(rtlim[1,]))
        col_leg <- colors[sel_raw]


        plot(
          0,
          xlab = "Time(min)",
          ylab = "Intensity",
          xlim = c(rtmin * 0.95 / 60, rtmax * 1.05 / 60),
          ylim = c(0, intmax * 1.05),
          type = "n",
          main = titles[i],
          cex.main = 0.9
        )

        ###We plot all the values
        for (j in seq_along(eics)) {
          if (is.null(eics[[j]][[i]]) || all(eics[[j]][[i]]$intensity == 0))
            next
          eic <- eics[[j]][[i]]
          lines(eic$time / 60,
                eic$intensity,
                lwd = 2,
                col = "black")
          lines(eic$time[eic$peak] / 60,
                eic$intensity[eic$peak],
                lwd = 2,
                col = colors[j])
        }
        legend(
          "topright",
          legend = names_raw[sel_raw],
          col = colors[sel_raw],
          lwd = 2,
          cex = 0.7
        )
      }
    }

  extracted_eics <-
    bpmapply(
      raw_files,
      peaktables,
      FUN = extract_EIC_raw_file,
      SIMPLIFY = FALSE,
      BPPARAM = bpp,
      MoreArgs = list(
        mz = vmz[found_signals],
        rt = vrt[found_signals],
        mz_margin = MARGIN_MZ,
        rt_margin = MARGIN_RT,
        num_detect = dm[, "num_detection"]
      )
    )

  extracted_RTs <- lapply(extracted_eics,"[[",i=2)
  extracted_eics <- lapply(extracted_eics,"[[",i=1)


  rts_val <- as.data.frame(do.call(cbind,extracted_RTs))
  colnames(rts_val) <- paste("rt","traw_files",sep="_")

  cnames <- colnames(rts_val)

  rts_val <- cbind(data.frame(name=titles[found_signals]),rts_val)
  colnames(rts_val) <- c("id",cnames)

  write.table(rts_val,file = OUTPUT_TARGETTED_RT_TABLE,row.names = FALSE,col.names = TRUE,sep = ";")

  ###we write the subDataMatrix
  write.table(dm[matched_features[found_signals],,drop=FALSE],file=OUTPUT_TARGETTED_INT_TABLE,sep=";",row.names = FALSE)
  found_signals <- !is.na(vmz)

  sink("/dev/null")
  pdf(OUTPUT_TARGET_PDF)
  plotPeaks(extracted_eics, titles = titles[found_signals], names_raw = traw_files)
  dev.off()
  sink(file = NULL)

}
