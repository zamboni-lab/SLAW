###Gap filling using an alignment model
suppressWarnings(suppressMessages(library(RSQLite, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(DBI, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(stringr, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))

###COnstant add columns name
##We add the isotopic pattern
ISO_DIST_COL <- "raw_isotopic_pattern"
ISO_NAME_COL <- "isotopic_pattern_annot"



##Argument passed by Python
args <- commandArgs(trailingOnly = TRUE)
#
# if(DEBUG){
# args <- c(
#   "U:/users/Alexis/data/slaw_evaluation/MTBLS1129/output_cluster/output_openms/processing_db.sqlite",
#   "U:/users/Alexis/data/datamatrix_31eb4257c5cdc54aabd59692e18f690c.csv",
#   "U:/users/Alexis/data/dm_output.csv",
#   "U:/users/Alexis/data/slaw_evaluation/MTBLS1129/output_cluster/output_openms/temp/alignement.rds",
#   "C:/Users/dalexis/Documents/dev/lcmsprocessing_docker/pylcmsprocessing/data/isotopes.tsv",
#   "4","3","intensity","0.001","15","0.01","5"
# )
# }


PATH_DB <- args[1]
PATH_DM <- args[2]
PATH_FILLED <- args[3]
PATH_MODEL <- args[4]
PATH_ISOTOPES <- args[5]
MAX_ISO <- as.integer(args[6])
MAX_CHARGE <- as.integer(args[7])
QUANT <- args[8]
MARGIN_MZ <- as.numeric(args[9])
PPM <- as.numeric(args[10])
DMZ <- as.numeric(args[11])
NUM_WORKERS <- as.integer(args[12])

TEMP_NAME <- paste(PATH_FILLED, 2, sep = "")

###This is jsut for evaluation
TEMP_FILLED <- paste(TEMP_NAME, "non_filled.csv", sep = "_")
nullvar <- file.copy(PATH_DM, TEMP_FILLED)

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




bpp <- NULL
if (get_os() == "win") {
  bpp <- SnowParam(workers = NUM_WORKERS, progressbar = TRUE)
} else{
  bpp <- MulticoreParam(workers = NUM_WORKERS, progressbar = TRUE)
}


dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
all_peaktables <-
  dbGetQuery(
    dbb,
    "SELECT output_ms FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1"
  )[, 1]
all_samples <-
  dbGetQuery(
    dbb,
    "SELECT path FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1"
  )[, 1]
# if(DEBUG){
# all_peaktables <- paste("U:/users/Alexis/data/slaw_evaluation/MTBLS1129/output_cluster/output_openms/OPENMS/peaktables/",basename(all_peaktables),sep="")
# all_samples <- paste("U:/users/Alexis/data/slaw_evaluation/MTBLS1129/mzML/",basename(all_samples),sep="")
# }
dbDisconnect(dbb)

isotopes_table <- fread(PATH_ISOTOPES, sep = " ")


cnames <- fread(PATH_DM,
                sep = ",",
                header = TRUE,
                nrows = 1)
cnames <- colnames(cnames)
dm <- fread(PATH_DM, sep = ",", select = c("mz", "rt"))
max_feat <- nrow(dm)
BY_BATCH <- 1000

rm(dm)
batches <- seq(1, max_feat, by = BY_BATCH)
batches[1] <- 0
if (batches[length(batches)] != max_feat) {
  batches <- c(batches, max_feat + 2)
} else{
  batches[length(batches)] <- max_feat + 2
}

for (idx in 1:(length(batches) - 1)) {
  dm <-
    fread(input = PATH_DM,
          nrows = batches[idx + 1] - batches[idx] - 1,
          skip = batches[idx])
  colnames(dm) <- cnames
  # dm <- fread(PATH_DM,sep=",")
  sel_columns <-
    c("min_mz", "max_mz", "mean_mz", "mean_rt", "mean_peakwidth")
  if (max(dm[["mean_peakwidth"]]) > max(dm[["mean_rt"]])) {
    dm[["mean_peakwidth"]] <- dm[["mean_peakwidth"]] / 60
  }
  dm_peaks <- dm[, ..sel_columns]
  ocnames <- colnames(dm)
  quant_prefix <-
    paste(str_split(ocnames[length(ocnames)], fixed("_"))[[1]][1], "_", sep =
            "")
  quant_cols <- which(startsWith(ocnames, quant_prefix))
  
  ###We reorder the data
  max_sample <- apply(dm[, ..quant_cols], 1, which.max)
  if (is.list(max_sample)) {
    max_sample <- sapply(max_sample, function(x) {
      x[1]
    })
  }
  by_sample <- split(seq_along(max_sample), max_sample)
  vord <- as.numeric(names(by_sample))
  isotopes_to_extract <-
    vector(mode = "list", length = length(quant_cols))
  for (ic in seq_along(isotopes_to_extract)) {
    isotopes_to_extract[[ic]] <- numeric(0)
  }
  isotopes_to_extract[vord] <- by_sample
  align <- readRDS(PATH_MODEL)
  rmz <- range(dm[["mean_mz"]])
  max_iso <- seq(rmz[1], rmz[2] + 200, by = 20)
  maxC <- ceiling(max_iso / 14)
  maxC <- maxC
  dist_iso <- sapply(maxC, function(x, miso, c13) {
    res <- rep(0, miso)
    for (is in 1:miso) {
      res[is] <- dbinom(is, x, c13)
      
    }
    return(res)
  }, miso = MAX_ISO, c13 = 0.01111)
  
  # message("Inferring",all_samples,"vs",,)
  ###We also extract the data to infer
  to_infer <-
    apply(dm[, ..quant_cols], 2, function(x) {
      which(is.na(x))
    })
  extractMissingInformations <-
    function(praw,
             peaks,
             isotopes,
             infer,
             align,
             table_iso,
             dist_c13,
             dm,
             quant,
             margin_mz = 0.001,
             max_iso = 4,
             max_charge = 2,
             ppm = 8,
             dmz = 0.002) {
      ###We create default vector to return if needed
      def_isotopes <- rep(NA, length(isotopes))
      def_missing <- rep(0.0, length(to_infer))
      
      suppressWarnings(suppressMessages(library(xcms, warn.conflicts = FALSE)))
      suppressWarnings(suppressMessages(library(pracma, warn.conflicts = FALSE)))
      suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))
      xraw <- suppressWarnings(suppressMessages(xcmsRaw(praw)))
      peaks <- fread(peaks)
      
      noise_level <- quantile(peaks$intensity, probs = 0.03)
      ort <- order(dm[infer, ][[4]])
      ###smal cheat the correction is the closest
      vmatch  <-
        xcms:::findEqualGreaterM(align[, 1], values = dm[infer[ort], ][[4]])
      align <- rbind(align, align[nrow(align), ])
      ###integrate the missing peak
      rt_corr <- dm[[4]][infer[ort]] + align[vmatch, 2]
      
      ###we first perform gap filling
      extractIntensity <-
        function(min_mz,
                 max_mz,
                 rt,
                 peakwidth,
                 idx,
                 quant,
                 xraw,
                 margin_mz) {
          rt_min <- rt - peakwidth / 2
          rt_max <- rt + peakwidth / 2
          min_mz <- min_mz - margin_mz
          max_mz <- max_mz + margin_mz
          tval <-
            tryCatch(
              rawEIC(
                xraw,
                rtrange = c(rt_min * 60, rt_max * 60),
                mzrange = c(min_mz, max_mz)
              ),
              error = function(e) {
                return(NA)
              }
            )
          if (length(tval) == 1 && is.na(tval))
            return(0)
          val <- 0
          if (quant == "intensity") {
            val <- trapz(x = xraw@scantime[tval[[1]]], tval[[2]])
          } else{
            val <- max(tval[[2]])
          }
          ##We extend the limit one more time.
          max_extension <- 5
          current_extension <- 0
          while (all(tval[[2]] == 0) &
                 current_extension < max_extension) {
            min_mz <- min_mz - margin_mz
            max_mz <- max_mz + margin_mz
            tval <-
              rawEIC(
                xraw,
                rtrange = c(rt_min * 60, rt_max * 60),
                mzrange = c(min_mz, max_mz)
              )
            if (any(tval[[2]] > 0)) {
              min_mz <- min_mz - 0.003
              max_mz <- max_mz + 0.003
              tval <-
                rawEIC(
                  xraw,
                  rtrange = c(rt_min * 60, rt_max * 60),
                  mzrange = c(min_mz, max_mz)
                )
            }
            current_extension <- current_extension + 1
          }
          
          if (quant == "intensity") {
            val <- trapz(x = xraw@scantime[tval[[1]]], tval[[2]])
          } else{
            val <- max(tval[[2]])
          }
          return(val)
        }
      inferred_values <-
        tryCatch(
          mapply(
            dm[["min_mz"]][infer[ort]],
            dm[["max_mz"]][infer[ort]],
            rt_corr,
            dm[["mean_peakwidth"]][infer[ort]],
            infer[ort],
            FUN = extractIntensity,
            MoreArgs = list(
              quant = quant,
              xraw = xraw,
              margin_mz = margin_mz
            )
          ),
          error = function(e) {
            return(def_missing)
          }
        )
      
      
      inferred_values <-
        inferred_values[order(ort, decreasing = FALSE)]
      
      ##We map the feature to the peaktable
      dm2 <- dm[isotopes, , drop = FALSE]
      visotopes <- list()
      if (length(isotopes) != 0) {
        mapFeatures <- function(dm2, pt) {
          mz_r <- max(dm2[, max_mz - min_mz])
          
          vmm <-
            xcms:::fastMatch(dm2[["mean_mz"]], pt[["mz"]], tol = mz_r)
          
          vnull <- which(sapply(vmm, is.null))
          while (length(vnull) > 0) {
            mz_r <- mz_r * 1.5
            vmm[vnull] <-
              xcms:::fastMatch(dm2[["mean_mz"]][vnull], pt[["mz"]], tol = mz_r)
            vnull <- which(sapply(vmm, is.null))
          }
          
          vv <- mapply(
            dm2[["mean_rt"]],
            vmm,
            FUN = function(x, y, rt_pp) {
              y[which.min(rt_pp[y])]
            },
            MoreArgs = list(rt_pp = pt[["rt"]])
          )
          return(vv)
        }
        
        vmap <- mapFeatures(dm2, peaks)
        
        ##Depending of the max isotopes we calculate the isotopic reation
        
        # if(FALSE){
        #   library(xcms)
        #
        #   dmm <- fread("U:/users/Alexis/data/slaw_evaluation/all_2mins/output/datamatrices/datamatrix_4a58048c5bb6c72f8a7dd470860d28cc.csv",sep=",")
        #   which.min(abs(dmm$mz-765.526968053277))
        #
        #   peaks <- "U:/users/Alexis/data/slaw_evaluation/all_2mins/output/OPENMS/peaktables/MSplate8_WT_NA_A12_plateY_rep1_SALT_GEII_1-09_792.csv"
        #   peaks <- fread(peaks,sep=",")
        #   ppeak <- peaks[2330,]
        #   peak <- as.numeric(ppeak)
        #   names(peak) <- names(ppeak)
        #   praw <- "U:/users/Alexis/data/slaw_evaluation/all_2mins/mzML/MSplate8_WT_NA_A12_plateY_rep1_SALT_GEII_1-09_792.mzML"
        #   xraw <- xcmsRaw(praw)
        #   iso_table <- read.table("C:/Users/dalexis/Documents/dev/lcmsprocessing_docker/pylcmsprocessing/data/isotopes.tsv",sep=" ",header=TRUE,stringsAsFactors = FALSE)
        #
        # }
        #
        extractIsotopes <- function(peak,
                                    xraw,
                                    iso_table,
                                    dist_iso,
                                    noise_level,
                                    max_iso = 4,
                                    max_charge = 2,
                                    ppm = 10,
                                    dmz = 0.004) {
          # print(iso_table)
          mz_peak <- peak["mz"]
          int_peak <- peak["height"]
          ##We select the closest feature
          sel_scan <-
            which.min(abs(peak["rt"] * 60 - xraw@scantime))
          mzrange <-
            c(peak[["mz_min"]] - 0.001, peak["mz_max"] + max(4, max_iso) + 0.5)
          sc <- getScan(xraw, sel_scan, mzrange = mzrange)
          iso_table[2:nrow(iso_table), "proportion"] <-
            iso_table[2:nrow(iso_table), "proportion"] * iso_table[2:nrow(iso_table), "max"]
          maxCp <- (peak["mz_max"] / 14)
          massC13 <- iso_table[["massdiff"]][1]
          massdiffC13 <- massC13 * (1:max_iso)
          idp <- floor(maxCp / 20)
          distC13 <- dist_iso[, idp]
          tol <- min((ppm * mz_peak * 1e-6), dmz)
          max_theo_int <- int_peak * distC13
          if (max(max_theo_int) < noise_level)
            return(NA)
          charge <- NA
          sel_iso <- NULL
          mass_diff <- NULL
          for (charge in max_charge:1) {
            mz_theo <- mz_peak + massdiffC13 / charge
            miso <- xcms:::fastMatch(mz_theo, sc[, 1], tol = tol)
            
            ###We filter by  intensity
            if (is.null(miso[[1]]))
              next
            ###We unroll the index of isotopic pattern
            idx <- 1
            while (idx <= max_iso &&
                   !is.null(miso[[idx]]) &&
                   sc[miso[[idx]][1], 2] <= max_theo_int[idx]) {
              miso[[idx]] <- miso[[idx]][1]
              idx <- idx + 1
            }
            if (idx == 1)
              return(NA)
            sel_iso <- sc[unlist(miso[1:(idx - 1)]), 2]
            mass_diff <- sc[unlist(miso[1:(idx - 1)]), 1] - mz_peak
            break
          }
          if (is.null(sel_iso))
            return(NA)
          sel_iso_c13 <- c(sc[1, 2], sel_iso)
          name_c13 <-
            c("", paste(1:length(sel_iso), rep("C13", length(sel_iso)), sep = ""))
          mass_diff <- c(0, mass_diff)
          df_to_return <-
            data.frame(
              name = name_c13,
              int = sel_iso_c13,
              massdiff = mass_diff,
              charge = rep(charge, length(mass_diff))
            )
          
          ##We remove the peaks which have already been matched to C13
          sel_peaks <- 1:nrow(sc)
          sel_peaks <- sel_peaks[-unlist(miso[1:(idx - 1)])]
          
          
          mz_other_isotopes <-
            iso_table[2:nrow(iso_table), "massdiff"]
          name_other_isotopes <-
            iso_table[2:nrow(iso_table), "isotope"]
          dist_other_isotopes <-
            iso_table[2:nrow(iso_table), "proportion"]
          supp_mz <- unlist(mz_other_isotopes / charge + mz_peak)
          msupp <-
            xcms:::fastMatch(supp_mz, sc[sel_peaks, 1], tol = tol, symmetric = FALSE)
          
          supp_found <- which(!sapply(msupp, is.null))
          if (length(supp_found) == 0) {
            return(df_to_return)
          } else{
            msupp[supp_found] <- lapply(msupp[supp_found], "[", i = 1)
            msupp[supp_found] <-
              lapply(msupp[supp_found], function(x, idx) {
                idx[x]
              }, idx = sel_peaks)
            ###We check intensity
            supp_found <-
              supp_found[sc[unlist(msupp[supp_found]), 2] < 1.1 * int_peak * dist_other_isotopes[supp_found]]
            if (length(supp_found) == 0)
              return(df_to_return)
            ###We check if there is a conflict if so we group it
            matched_features <- msupp[supp_found]
            supp_table <-
              data.frame(
                name = unlist(name_other_isotopes[supp_found]),
                int = sc[unlist(msupp[supp_found]), 2],
                massdiff = sc[unlist(msupp[supp_found]), 1] -
                  mz_peak,
                charge = rep(charge, length(supp_found))
              )
            supp_table <-
              by(
                supp_table,
                INDICES = supp_table$int,
                FUN = function(x) {
                  data.frame(
                    name = paste(x[["name"]], collapse = "@"),
                    int = x[["int"]][1],
                    massdiff = x[["massdiff"]][1],
                    charge = x[["charge"]][1]
                  )
                }
              )
            supp_table <- do.call(rbind, supp_table)
            
            if (length(supp_found) == 0)
              return(df_to_return)
            df_to_return <- rbind(df_to_return, supp_table)
            return(df_to_return)
          }
        }
        visotopes <-
          tryCatch(
            apply(
              peaks[vmap],
              1,
              extractIsotopes,
              xraw = xraw,
              iso_table = table_iso,
              dist_iso = dist_c13,
              noise_level = noise_level,
              max_iso = max_iso,
              max_charge = max_charge,
              ppm = ppm,
              dmz = dmz
            ),
            error = function(e) {
              return(def_isotopes)
            }
          )
      }
      return(list(inferred_values, visotopes))
    }
  #
  # vmap <- vector(mode="list",length=length(all_samples))
  # if(DEBUG){
  # for(aidx in seq_along(all_samples)){
  #   print(aidx)
  # vmap[[aidx]] <- extractMissingInformations(all_samples[[aidx]],
  #                            all_peaktables[[aidx]],
  #                            isotopes_to_extract[[aidx]],
  #                            to_infer[[aidx]],
  #                            align@rt_correction[[aidx]],
  #                            dm=dm_peaks,quant=QUANT,table_iso=isotopes_table,
  # dist_c13=dist_iso,margin_mz=MARGIN_MZ,max_iso = MAX_ISO,
  # max_charge = MAX_CHARGE, ppm = PPM,dmz = DMZ)
  # }
  
  # if(!DEBUG){
  vmap <-
    bpmapply(
      all_samples,
      all_peaktables,
      isotopes_to_extract,
      to_infer,
      align@rt_correction,
      FUN = extractMissingInformations,
      MoreArgs = list(
        dm = dm_peaks,
        quant = QUANT,
        table_iso = isotopes_table,
        dist_c13 =
          dist_iso,
        margin_mz = MARGIN_MZ,
        max_iso = MAX_ISO,
        max_charge = MAX_CHARGE,
        ppm = PPM,
        dmz = DMZ
      ),
      SIMPLIFY = FALSE,
      BPPARAM = bpp
    )
  # }
  ###We fill the column
  for (iquant in seq_along(quant_cols)) {
    if (length(to_infer[[iquant]]) == 0)
      next
    sel_col <- quant_cols[iquant]
    dm[[sel_col]][to_infer[[iquant]]] <- vmap[[iquant]][[1]]
  }
  
  
  name_col <- rep(NA_character_, nrow(dm))
  dist_col <- rep(NA_real_, nrow(dm))
  
  for (is in seq_along(vmap)) {
    if (length(isotopes_to_extract[[is]]) == 0)
      next
    names <- sapply(vmap[[is]][[2]], function(x) {
      if (!is.data.frame(x))
        return(NA)
      sel_idx <- 2:nrow(x)
      n0 <- paste("M")
      n_supp <-
        paste("M+",
              sprintf("%0.4f", x[sel_idx, "massdiff"]),
              "_",
              x[sel_idx, "name"],
              sep = "",
              collapse = "|")
      return(paste(n0, n_supp, sep = "|"))
    })
    dists <- sapply(vmap[[is]][[2]], function(x) {
      if (!is.data.frame(x))
        return(NA)
      paste(sprintf("%0.4f", x[, "int"] / max(x[, "int"])), collapse = "|")
    })
    
    ###DEBUG ONLY
    ###
    # is_buggued <- sapply(vmap[[is]][[2]],function(x){
    #   if(!is.data.frame(x)) return(FALSE)
    #   return(any(x[1:(nrow(x)-1),"int"]<x[2:nrow(x),"int"]))
    # })
    #
    # if(any(is_buggued)){
    #   for(pp in which(is_buggued)){
    #   message("pos:",isotopes_to_extract[[is]][pp],"|infos:",dm[isotopes_to_extract[[is]][pp],1:2],
    #           "|file:",all_samples[is],all_peaktables[is])
    #   }
    # }
    
    name_col[isotopes_to_extract[[is]]] <- names
    dist_col[isotopes_to_extract[[is]]] <- dists
  }
  ###The data matrix is expanded
  infos_idx <- 1:(quant_cols[1] - 1)
  quant_idx <- quant_cols
  new_names <-
    c(colnames(dm)[infos_idx], ISO_NAME_COL, ISO_DIST_COL, colnames(dm)[quant_idx])
  quant_cols_2 <- quant_cols + 2
  dm <-
    cbind(dm[, ..infos_idx], name_col, dist_col, dm[, ..quant_idx])
  colnames(dm) <- new_names
  ww <- fwrite(dm, TEMP_NAME, append = TRUE)
}
message("Gap filling finished")
ww <- file.rename(PATH_DM, PATH_FILLED)
ww <- file.rename(TEMP_NAME, PATH_DM)
