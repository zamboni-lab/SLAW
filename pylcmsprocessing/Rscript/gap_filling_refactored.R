###Gap filling using an alignment model
suppressWarnings(suppressMessages(library(RSQLite, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(DBI, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(stringr, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(BiocParallel, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(rhdf5, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(pracma, warn.conflicts = FALSE)))
suppressWarnings(suppressMessages(library(xcms, warn.conflicts = FALSE)))

### DEBUG
DEBUG <- FALSE

###COnstant add columns name
##We add the isotopic pattern
ISO_DIST_COL <- "isotopic_pattern_rel"
ISOABS_DIST_COL <- "isotopic_pattern_abs"
ISO_NAME_COL <- "isotopic_pattern_annot"
EMPTY_ISOTOPES <- "EMPTY"
MAX_FILES <- 10

sink(file=stdout())
##Argument passed by Python
args <- commandArgs(trailingOnly = TRUE)

# testing
if (length(args)<10) {
  DEBUG <- TRUE
  DEBUG_OUTPUT <- "D:/Data/test_pwiz3023_slaw/"
  DEBUG_INPUT <- "D:/Data/test_pwiz3023/"
  args <- c("/output/processing_db.sqlite",
            "/output/temp/data_prefill_9a08a61020295d2c79a91ecec29a4892.csv",
            "/output/temp","/output/temp/alignment.rds","/output/temp/gap_filling.hdf5",
            "D:\\SW\\SLAW\\pylcmsprocessing\\data\\isotopes.tsv","quant"
            ,8,3,15.0,0.01,5,0.5,16)
  args <- sapply(args,str_replace,"/output/",DEBUG_OUTPUT)
}

PATH_DB <- args[1]
PATH_DM <- gsub('data_filled_','data_prefill_',args[2]) #str_replace(args[2],'data_filled_','data_prefill_')
PATH_TEMP <- args[3]
PATH_MODEL <- args[4]
HDF5_FILE <- args[5]
PATH_ISOTOPES <- args[6]
QUANT <- 'quant' #args[7]
MAX_ISO <- as.integer(args[8])
MAX_CHARGE <- as.integer(args[9])
PPM <- as.numeric(args[10])
DMZ <- as.numeric(args[11])
NUM_FILES <- as.integer(args[12])
FRAC_CUTOFF <- as.numeric(args[13])
NUM_WORKERS <- as.integer(args[13])

if(file.exists(gsub('data_prefill_','data_filled_',PATH_DM))){
  print("Matrix already exists. Skipping gap filling. Delete the file to repeat this step.") 
  # stop execution
  quit(save = "no", status = 0, runLast = FALSE)  
}
#Path of the temporary filled matrix
# TEMP_FILLED <- file.path(PATH_TEMP,"filled_dm.csv")
# TEMP_TRANSFER <- file.path(PATH_TEMP,"transfer_dm.csv")

#Path of missing informations
HDF5_FILE <- file.path(PATH_TEMP,"missing_infos.hdf5")

#Isotopes table
isotopes_table <- fread(PATH_ISOTOPES, sep = " ")

#Column names of the data matrix
cnames <- fread(PATH_DM,
                sep = "\t",
                header = TRUE,
                nrows = 1)
cnames <- colnames(cnames)
dm <- fread(PATH_DM, sep = "\t")
num_features <- nrow(dm)

#Handling of paralleliszation depending of OS
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

FEATURES_GROUP <- "FEATURES"
ISOTOPES_GROUP <- "ISOTOPES"
MISSING_FEATURES_PREFIX <- "MFEAT"
ISOTOPES_TO_FIND_PREFIX <- "FISOS"
ISOTOPES_TO_FIND_VALUES_PREFIX <- "FISOS_VALUES"
MISSING_FEATURES_VALUES_PREFIX <- "MFEAT_VALUES"


#We extract the informations on the samples
dbb <- dbConnect(RSQLite:::SQLite(), PATH_DB)
all_infos <-
  dbGetQuery(
    dbb,
    "SELECT path,output_ms,types FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1"
  )
all_peaktables <- all_infos[,2]
all_samples <- all_infos[,1]
if (DEBUG) {
  all_peaktables <- sapply(all_peaktables,str_replace,"/output/",DEBUG_OUTPUT)
  all_samples <- sapply(all_samples,str_replace,"/input/",DEBUG_INPUT)
}

#The number of cores is potentially readjusted if the raw files are too big
num_samples <- length(all_samples)

# Extracting quantitative columns and infos from the data matrix
quant_prefix <-
  paste(str_split(cnames[length(cnames)], fixed("_"))[[1]][1], "_", sep = "")
quant_cols <- which(startsWith(cnames, quant_prefix))
quant_mean <- apply(dm[,..quant_cols],1,mean,na.rm=TRUE)

# count the frequency of missing values in each row of dm
quant_na_frequency <- apply(dm[,..quant_cols], 1, function(x) sum(is.na(x)))/length(quant_cols)
features_to_fill <- which(quant_na_frequency >= FRAC_CUTOFF)

#We define the batches of files if thee is too much
seq_samples <- seq(0,num_features*length(quant_cols),by=700000)
seq_samples <- seq_samples%/%num_features
seq_samples[1] <- 0
if (seq_samples[length(seq_samples)] != num_samples) {
  seq_samples <- c(seq_samples, num_samples)
}

#We remove the storage if it exists
if(file.exists(HDF5_FILE)){
  dummy <-file.remove(HDF5_FILE)
}
dummy <- h5createFile(HDF5_FILE)
dummy <- h5createGroup(HDF5_FILE,FEATURES_GROUP)
dummy <- h5createGroup(HDF5_FILE,ISOTOPES_GROUP)

#We first check the id by batch
max_feat <- nrow(dm)
BY_BATCH <- 10000

batches <- seq(1, max_feat, by = BY_BATCH)
if (batches[length(batches)] != max_feat) {
  batches <- c(batches, max_feat)
} else{
  batches[length(batches)] <- max_feat + 1
}

#Extracting the features to infer from each samples
raw_isotopes_samples <- rep(0,num_features)
raw_isotopes_intensities <- rep(0,num_features)

for(isample in 1:(length(seq_samples)-1)){
  batch_idx <- (seq_samples[isample]+1):seq_samples[isample+1]
  dm_batch_idx <- quant_cols[batch_idx]
  missing <- apply(dm[,..dm_batch_idx],2,function(x){which(is.na(x))})
  if(is.matrix(missing)){
    missing <- split(missing,1:ncol(missing))
  }
  #We insert them into the hdf5
  for(idx_sample in seq_along(missing)){
    #We also write them by batch
    for(idx_batch in 1:(length(batches)-1)){#This is the length of the batches -1
      current_missing <- missing[[idx_sample]]
      # added to skip filling of rare features
      current_missing <- intersect(current_missing, features_to_fill) 
      sel_pos <- which(batches[idx_batch]<=current_missing&batches[idx_batch+1]>current_missing)
      feat_prefix <- paste(FEATURES_GROUP,"/",MISSING_FEATURES_PREFIX,batch_idx[idx_sample],"_",idx_batch,sep="")
      h5write(current_missing[sel_pos]-batches[idx_batch]+1,HDF5_FILE,feat_prefix)
    }
  }
  #We get the potential isotopes
  max_isos <- apply(dm[,..dm_batch_idx],1,function(x){
    pmax <- suppressWarnings(which.max(x))
    #Feature not found in this batch
    if(length(pmax)==0){
      return(c(0,0))
    }else{
      return(c(x[pmax],pmax))
    }
  })

  found_isotopes <- which(max_isos[1,] != 0)
  if(length(found_isotopes)>0){
    max_isos_val <- max_isos[1,found_isotopes]
    max_isos_sample <- batch_idx[max_isos[2,found_isotopes]]

    #We update the position vector if needed
    to_change <- which(raw_isotopes_intensities[found_isotopes]<max_isos_val)
    if(length(to_change)>0){
      raw_isotopes_intensities[found_isotopes[to_change]] <- max_isos_val[to_change]
      raw_isotopes_samples[found_isotopes[to_change]] <- max_isos_sample[to_change]
    }
  }
}

#We write the isotopic ratio in the hdf5 for each sample
for(samp in 1:num_samples){
  sel_isos <- which(raw_isotopes_samples==samp)

  for(idx_batch in 1:(length(batches)-1)){
    sel_pos <- which(batches[idx_batch]<=sel_isos&batches[idx_batch+1]>sel_isos)
    iso_prefix <- paste(ISOTOPES_GROUP,"/",ISOTOPES_TO_FIND_PREFIX,samp,"_",idx_batch,sep="")
    h5write(sel_isos[sel_pos]-batches[idx_batch]+1,HDF5_FILE,iso_prefix)
  }
}

mapFeatures <- function(dm2, pt) {
  mz_r <- min(max(dm2[, mz_max - mz_min]),0.03)
  
  vmm <- xcms:::fastMatch(dm2[["mz_mean"]], pt[["mz"]], tol = mz_r)
  
  vnull <- which(sapply(vmm, is.null))
  lim <- mz_r*5
  while (length(vnull) > 0 & mz_r<lim) {
    mz_r <- mz_r * 1.5
    vmm[vnull] <-
      xcms:::fastMatch(dm2[["mz_mean"]][vnull], pt[["mz"]], tol = mz_r)
    vnull <- which(sapply(vmm, is.null))
  }
  vfound <- which(!sapply(vmm, is.null))
  vv <- rep(NA,length(vmm))
  
  vv[vfound] <- sapply(
    vmm[vfound],
    FUN = function( y, rt_pp) {
      y[[which.min(rt_pp[y])]]
    },
    rt_pp = pt[["rt"]]
  )
  return(vv)
}

extractIntensity <-
  function(min_mz,
           max_mz,
           rt,
           peakwidth,
           idx,
           expected_intensity,
           quant,
           xraw,
           margin_mz,
           margin_ppm,
           margin_dmz) {
    
    integ_fun <- function(x, y) {
      trapz(x, y)
    }
    
    if (quant != "intensity") {
      
      integ_fun <- function(x, y) {
        max(y)
      }
      
    }
    
    rt_min <- rt - peakwidth / 2
    rt_max <- rt + peakwidth / 2
    #We shift the peak if it is negativer
    if (rt_min < 0) {
      rt_max <- rt_max + abs(rt_min)
      rt_min <- 0
    }
    
    
    mean_mz <- (min_mz + max_mz) / 2
    
    ###Margin in ppm.
    margin_mz <- max(margin_dmz, mean_mz * margin_ppm * 1e-6)
    min_mz <- min_mz - margin_mz
    max_mz <- max_mz + margin_mz
    total_extension <- margin_mz
    
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
      return(NA) # it was zero
    val <- NA
    ##We extend the limit one more time.
    max_extension <- 0.07
    current_extension <- 0
    old_int <- -1
    current_int <- integ_fun(x = xraw@scantime[tval[[1]]], tval[[2]])
    ###We keep expanding the trace mass.
    while ((abs(expected_intensity - current_int) <= abs(expected_intensity - old_int)) &
           (old_int != current_int) &
           (current_int < expected_intensity) &
           (total_extension < max_extension)) {
      min_mz <- min_mz - margin_mz
      max_mz <- max_mz + margin_mz
      tval <-
        rawEIC(
          xraw,
          rtrange = c(rt_min * 60, rt_max * 60),
          mzrange = c(min_mz, max_mz)
        )
      old_int <- current_int
      current_int <- integ_fun(x = xraw@scantime[tval[[1]]], tval[[2]])
      total_extension <- total_extension + margin_mz
    }
    val <- old_int
    if (val == -1) val <- current_int
    return(val)
  }

extractIsotopes <- function(peak,
                            xraw,
                            iso_table,
                            dist_iso,
                            noise_level,
                            max_iso = 4,
                            max_charge = 2,
                            ppm = 10,
                            dmz = 0.004) {
  mz_peak <- peak["mz"]
  int_peak <- peak["height"]
  ##We select the closest feature
  sel_scan <- which.min(abs(peak["rt"] * 60 - xraw@scantime))
  mzrange <- c(peak[["mz_min"]] - 0.001, peak["mz_max"] + max(4, max_iso) + 0.5)
  sc <- getScan(xraw, sel_scan, mzrange = mzrange)
  iso_table[2:nrow(iso_table), "proportion"] <- iso_table[2:nrow(iso_table), "proportion"] * iso_table[2:nrow(iso_table), "max"]
  maxCp <- (peak["mz_max"] / 20)
  massC13 <- iso_table[["massdiff"]][1]
  massdiffC13 <- massC13 * (1:max_iso)
  # idp <- floor(maxCp) # /20
  # if (idp == 0) {
  #   return(NA)
  # }
  # distC13 <- dist_iso[, idp]
  # max_theo_int <- int_peak * distC13
  # if (max(max_theo_int) < noise_level) {
  #   return(NA)
  # }
  # 
  charge <- NA
  sel_iso <- NULL
  mass_diff <- NULL
  tol <- min((ppm * mz_peak * 1e-6), dmz)
  for (charge in max_charge:1) {
    mz_theo <- mz_peak + massdiffC13 / charge
    miso <- xcms:::fastMatch(mz_theo, sc[, 1], tol = tol)
    
    ###We filter by  intensity
    if (is.null(miso[[1]]))
      next
    
    ###We unroll the index of isotopic pattern
    idx <- 1
    while (idx <= max_iso &&
           !is.null(miso[[idx]])) {
      miso[[idx]] <- miso[[idx]][1]
      idx <- idx + 1
    }
    # DEBUG: REMOVE FILTER ON MAX_THEO_INT of single isotopes >> doesn't work for heavy stuff.
    # while (idx <= max_iso &&
    #        !is.null(miso[[idx]]) &&
    #        sc[miso[[idx]][1], 2] <= max_theo_int[idx]) {
    #   miso[[idx]] <- miso[[idx]][1]
    #   idx <- idx + 1
    # }
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
  } else {
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

extractMissingInformation <- function(praw,
                                      peaks,
                                      isotopes,
                                      infer,
                                      align,
                                      table_iso,
                                      dist_c13,
                                      dm,
                                      quant,
                                      quant_mean,
                                      hdf5_file,
                                      margin_dmz = 0.003,
                                      margin_ppm = 10,
                                      margin_mz = 0.001,
                                      max_iso = 4,
                                      max_charge = 2,
                                      ppm = 8,
                                      dmz = 0.002,
                                      detect_isotopes = TRUE) {
  
  ###We create default vector to return if needed
  def_isotopes <- rep(NA, length(isotopes))
  def_missing <- rep(NA, length(infer))
  
  #Isotopes and inference
  xraw <- praw
  if (is.character(praw)) {
    xraw <- suppressWarnings(suppressMessages(xcmsRaw(praw)))
  }
  if (is.character(peaks)) {
    peaks <- fread(peaks)
  }
  noise_level <- 1000 # quantile(peaks$intensity, probs = 0.005) ### THE CUTOFF HAS MAJOR IMPLICATIONS. ABSOLUTE OR QUANTILE? CHANGED FROM 3% to 0.5%
  noise_level_height <- 10 # 1000 quantile(peaks$height, probs = 0.005)
  ###If the two values are equl  the noise is read from a file
  if (noise_level_height == noise_level) {
    vsc <- getScan(xraw, floor(length(xraw@scantime) / 2))[, 2]
    noise_level_height <- quantile(vsc, probs = 0.005)
  }
  
  ort <- order(dm[infer,][[4]])
  ###small cheat the correction is the closest
  all_corrs <- rep(0, length(infer[ort]))
  if (is.na(align[1, 1])) {
    dev_rt <- rep(0, length(infer[ort]))
  }else {
    vmatch <-
      xcms:::findEqualGreaterM(align[, 1], values = dm[infer[ort],][[4]])
    align <- rbind(align, align[nrow(align),])
    dev_rt <- align[vmatch, 2]
  }
  ###integrate the missing peak
  rt_corr <- dm[[4]][infer[ort]] + dev_rt
  
  ###we first perform gap filling
  inferred_values <-
    tryCatch(
      mapply(
        dm[["mz_min"]][infer[ort]],
        dm[["mz_max"]][infer[ort]],
        rt_corr,
        dm[["peakwidth_mean"]][infer[ort]],
        infer[ort],
        quant_mean[infer[ort]],
        FUN = extractIntensity,
        MoreArgs = list(
          quant = quant,
          xraw = xraw,
          margin_mz = margin_mz,
          margin_ppm = margin_ppm,
          margin_dmz = margin_dmz
        )
      ),
      error = function(e) {
        print(attr(e, "traceback"))
        print(e)
      }
    )
  inferred_values <- inferred_values[order(ort, decreasing = FALSE)]
  
  visotopes_pos <- numeric(0)
  visotopes_values_pos <- list()
  visotopes_idx <- numeric(0)
  if (detect_isotopes & length(isotopes) != 0) {
    # if (24 %in% isotopes)
    # {
    #   print(xraw@filepath[1])
    # }
    
    ##We map the feature to the peaktable
    dm2 <- dm[isotopes, , drop = FALSE]
    visotopes <- list()
    if (length(isotopes) != 0) {
      
      vmap <- mapFeatures(dm2, peaks)
      peaks_to_find <- which(!is.na(vmap))
      
      if (length(peaks_to_find) != 0) {
        #Incorrect mapping case.
        visotopes <-
          tryCatch(
            apply(
              peaks[vmap[peaks_to_find]],
              1,
              extractIsotopes,
              xraw = xraw,
              iso_table = table_iso,
              dist_iso = dist_c13,
              noise_level = noise_level_height,
              max_iso = max_iso,
              max_charge = max_charge,
              ppm = ppm,
              dmz = dmz),
            error = function(e) {
              print(e)
              # print("peaks_to_find")
              # print(peaks_to_find)
              # print("vmap")
              # print(vmap)
              # print("vmap[peaks_to_find]")
              # print(vmap[peaks_to_find])
            }
          )
      }
    }
    
    if (length(visotopes) != 0) { #Can be 0 if the mapping fail
      correct_pos <- which(sapply(visotopes, is.data.frame))
      visotopes_pos <- peaks_to_find[correct_pos]
      visotopes_values_pos <- visotopes[correct_pos]
      visotopes_idx <- isotopes[visotopes_pos]
    }
    
  }
  #We only store the inferred values to avoid the cost
  inferred_values_pos <- which(inferred_values != 0.0)
  ##Both inferred value will be written in the data
  return(list(list(visotopes_pos, visotopes_values_pos, visotopes_idx),
              list(inferred_values_pos, inferred_values[inferred_values_pos])))
}

extractMissingInformationHDF5 <- function(praw,
                                           peaks,
                                           path_isotopes, #hdf5 file path
                                           path_infer, #hdf5 file path
                                           align,
                                           table_iso,
                                           dist_c13,
                                           dm,
                                           quant,
                                           quant_mean,
                                           hdf5_file,
                                           extractMissingInformationFun,
                                           margin_dmz = 0.003,
                                           margin_ppm = 10,
                                           margin_mz = 0.001,
                                           max_iso = 4,
                                           max_charge = 2,
                                           ppm = 8,
                                           dmz = 0.002,
                                           detect_isotopes = TRUE) {
  
  isotopes <- h5read(hdf5_file, path_isotopes)
  infer <- h5read(hdf5_file, path_infer)
  rlist <- extractMissingInformationFun(praw = praw,
                                        peaks = peaks,
                                        isotopes,
                                        infer,
                                        align = align,
                                        table_iso = table_iso,
                                        dist_c13 = dist_c13,
                                        dm = dm,
                                        quant = quant,
                                        quant_mean = quant_mean,
                                        hdf5_file = hdf5_file,
                                        margin_dmz = margin_dmz,
                                        margin_ppm = margin_ppm,
                                        margin_mz = margin_mz,
                                        max_iso = max_iso,
                                        max_charge = max_charge,
                                        ppm = ppm,
                                        dmz = dmz,
                                        detect_isotopes = detect_isotopes)

  return(rlist)
}

extractMissingInformationAllBatches <- function(praw,
                                                 peaks,
                                                 paths_isotopes, #hdf5 file path
                                                 paths_infer, #hdf5 file path
                                                 align,
                                                 extractMissingInformationFun, ##FOr parallelization purpose
                                                 table_iso,
                                                 dist_c13,
                                                 dm,
                                                 quant,
                                                 quant_mean,
                                                 hdf5_file,
                                                 margin_dmz = 0.003,
                                                 margin_ppm = 10,
                                                 margin_mz = 0.001,
                                                 max_iso = 4,
                                                 max_charge = 2,
                                                 ppm = 8,
                                                 dmz = 0.002,
                                                 detect_isotopes = TRUE) {


  ##Same function but in this case we read the file once and just pass it around across the dataset
  praw <- suppressWarnings(suppressMessages(xcmsRaw(praw)))
  peaks <- fread(peaks)
  all_outputs <- mapply(paths_isotopes, paths_infer, FUN = extractMissingInformationHDF5,
                        MoreArgs = list(praw = praw,
                                        peaks = peaks,
                                        align = align,
                                        table_iso = table_iso,
                                        dist_c13 = dist_c13,
                                        extractMissingInformationFun = extractMissingInformation,
                                        dm = dm,
                                        quant = quant,
                                        quant_mean = quant_mean,
                                        hdf5_file = hdf5_file,
                                        margin_dmz = margin_dmz,
                                        margin_ppm = margin_ppm,
                                        margin_mz = margin_mz,
                                        max_iso = max_iso,
                                        max_charge = max_charge,
                                        ppm = ppm,
                                        dmz = dmz,
                                        detect_isotopes = detect_isotopes), SIMPLIFY = FALSE)

  res_isos <- lapply(all_outputs, "[[", i = 1)
  res_infers <- lapply(all_outputs, "[[", i = 2)
  return(list(res_isos, res_infers))
}


#This part select the files and perform the optimization if needed.
optim_idx <- seq_along(all_samples)

if("QC" %in% all_infos[,3]){
  optim_idx <- which(str_detect(all_infos[,3],fixed("QC")))
}

if(length(optim_idx)>NUM_FILES){
  sel_idx <- sample(seq_along(optim_idx),size=NUM_FILES)
  optim_idx <- optim_idx[sel_idx]
}

dbDisconnect(dbb)

isotopes_table <- fread(PATH_ISOTOPES, sep = " ")
cnames <- fread(PATH_DM,
                sep = "\t",
                header = TRUE,
                nrows = 1)
cnames <- colnames(cnames)
colnames(dm) <- cnames

sel_columns <-
  c("mz_min", "mz_max", "mz_mean", "rt_mean", "peakwidth_mean")
if (max(dm[["peakwidth_mean"]]) > max(dm[["rt_mean"]])) {
  dm[["peakwidth_mean"]] <- dm[["peakwidth_mean"]] / 60
}
dm_peaks <- dm[, ..sel_columns]

align <- readRDS(PATH_MODEL)
quant_cols_optim <- quant_cols[optim_idx]

# Tfis is just a small subset on which to optimize the gap filling
isotopes_to_extract <-
  vector(mode = "list", length = length(quant_cols))
for (ic in seq_along(isotopes_to_extract)) {
  isotopes_to_extract[[ic]] <- numeric(0)
}
dist_iso <- NULL

# We also extract the data to infer
to_infer <-
  apply(dm[, ..quant_cols_optim], 2, function(x,max_infer=1000) {
    to_infer <- which(is.na(x))
    if(length(to_infer)>max_infer){
      to_infer <- sample(to_infer,max_infer)
    }
    to_infer
  })

if(is.matrix(to_infer)){
  to_infer <- split(to_infer,1:ncol(to_infer))
}

optim_intensity <- apply(dm[,..quant_cols_optim],1,mean,na.rm=TRUE)




###We keep increasing the tolerance until the distribution stat to match.
optimizeParameters <-function(praws,peaks,isotopes,infer,
                             align,table_iso,dist_c13,
                             dm,quant,optim_intensity,
                             max_iso = 4,max_charge = 2,
                             ppm = 8,dmz = 0.002, by_file = 300,num_points=5){

  #We don t use the features withot intensity values
  to_use <- which(!is.nan(optim_intensity))

  ###We select at most 100 values by files.
  infer <- lapply(infer,FUN = function(x,max_size,to_keep){
    x <- x[x %in% to_keep]
    if(length(x)>max_size){
      return(sample(x,max_size))
    }
    return(x)
  },max_size=by_file,to_keep=to_use)

  ###Reading the file only once
  xraws <- suppressWarnings(suppressMessages(sapply(praws,xcmsRaw)))

  ###Onece the best margin has been found we can recompute the distribution
  score_gap_filling <- function(x,xraws,
                                peaks,
                                isotopes,
                                infers,
                                aligns,
                                table_iso,
                                dist_c13,
                                dm,
                                quant,
                                optim_intensity,
                                max_iso = 4,
                                max_charge = 2,
                                ppm = 8,
                                dmz = 0.002){

    vcomps <- mapply(xraws,peaks,isotopes,
                     infers,aligns, FUN = extractMissingInformation,
                     MoreArgs=list(table_iso=table_iso,dist_c13=dist_c13,dm=dm,quant=quant,quant_mean=optim_intensity,
                                   margin_mz = x,
                                   margin_ppm = x[2],
                                   margin_dmz = x[1],
                                   max_iso = max_iso,max_charge = max_charge,ppm = ppm,
                                   dmz = dmz,detect_isotopes = FALSE),SIMPLIFY = FALSE)

    score <- 0
    for(idx in seq_along(xraws)){
      ref_int <- optim_intensity[infers[[idx]]]
      ref_int <- ref_int[ vcomps[[idx]][[2]][[1]]]
      opt_int <- vcomps[[idx]][[2]][[2]]
      to_average <- abs(log10(opt_int+2)-log10(ref_int+2))
      score <- score + mean(to_average[optim_intensity[infers[[idx]]]!=0],na.rm = TRUE)
    }
    score <- score/length(xraws)
    return(score)
  }


  #This is potentially very costly in memory so we reduce it
  dmz_seq <- exp(seq(log(0.0005),log(0.05),length=num_points))
  ppm_seq <- exp(seq(log(2),log(30),length=num_points))
  grid_seq <- expand.grid(dmz_seq,ppm_seq)
  veval <- apply(grid_seq,1,FUN = score_gap_filling,xraws=xraws,
        peaks=peaks,isotopes=isotopes,infer=infer,
        align=align,table_iso=table_iso,dist_c13=dist_c13,
        dm=dm,quant=quant,optim_intensity=optim_intensity,max_iso = max_iso,
        max_charge = max_charge,ppm = ppm,
        dmz = dmz)

  best_par <- unlist(grid_seq[which.min(veval),])
  best_val <- min(veval)


  ###one we found a local minima, for the three minimal value we just do a gradient descent.
  dmz_step <- 0.0001
  ppm_step <- 0.2
  new_val <- -1
  new_par <- unlist(best_par)
  MAX_ITS <- 10
  current_it <- 1
  while(new_val<best_val & current_it < MAX_ITS){

    if(new_val<best_val & new_val!=-1){
      best_val <- new_val
      best_par <- new_par
    }

    val_dmz <- score_gap_filling(c(best_par[1]+dmz_step,best_par[2]),xraws=xraws,
    peaks=peaks,isotopes=isotopes,infer=infer,
    align=align,table_iso=table_iso,dist_c13=dist_c13,
    dm=dm,quant=quant,optim_intensity=optim_intensity,max_iso = max_iso,
    max_charge = max_charge,ppm = ppm,
    dmz = dmz)

    grad_dmz <- (val_dmz-best_val)/dmz_step

    val_ppm <- score_gap_filling(c(best_par[[1]],best_par[[2]]+ppm_step),xraws=xraws,
                                 peaks=peaks,isotopes=isotopes,infer=infer,
                                 align=align,table_iso=table_iso,dist_c13=dist_c13,
                                 dm=dm,quant=quant,optim_intensity=optim_intensity,max_iso = max_iso,
                                 max_charge = max_charge,ppm = ppm,
                                 dmz = dmz)
    grad_ppm <- (val_ppm-best_val)/ppm_step
    if(grad_ppm==0){
      new_ppm <- best_par[2]+runif(1,-1,1)*ppm_step
    }else{
      new_ppm <- best_par[2]-grad_ppm*ppm_step
    }

    if(grad_dmz==0){
      new_dmz <- best_par[2]+runif(1,-1,1)*dmz_step
    }else{
      new_dmz <- best_par[2]-grad_dmz*dmz_step
    }

    new_par <- c(new_dmz,new_ppm)
    new_val <- score_gap_filling(new_par,xraws=xraws,
                      peaks=peaks,isotopes=isotopes,infer=infer,
                      align=align,table_iso=table_iso,dist_c13=dist_c13,
                      dm=dm,quant=quant,optim_intensity=optim_intensity,max_iso = max_iso,
                      max_charge = max_charge,ppm = ppm,
                      dmz = dmz)
    current_it <- current_it+1
  }
  return(c(best_par,best_val))
}

##We only keep the best parameters
temp_par <- tryCatch(optimizeParameters(all_samples[optim_idx],all_peaktables[optim_idx],isotopes_to_extract[optim_idx],
                               to_infer,align@rt_correction[optim_idx],
                               table_iso = isotopes_table,
                               dist_c13 = dist_iso,
                               dm = dm_peaks, quant = QUANT,
                               optim_intensity = optim_intensity,
                               max_iso = 4, max_charge = 2,
                               ppm = 8, dmz = 0.002, num_points = 2),
                               error=function(e){print(e);print(attr(e, "traceback"));return(c(0.005,10))})

margin_dmz <- temp_par[1]
margin_ppm <- temp_par[2]
rm(dm)

#We read the elements that we will use for gap-filling
rmz <- dm_peaks[["mz_mean"]][c(1,length(dm_peaks[["mz_mean"]]))]
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

# In this step the files are opened sequentially.
for(isamp in seq_along(all_peaktables)){
  #Extract all the missing values for this particular samples

  #The prefix to READ the infos from the hdf5
  feat_input_infos <- paste(FEATURES_GROUP,"/",MISSING_FEATURES_PREFIX,isamp,"_",1:(length(batches) - 1),sep="")
  isos_input_infos <- paste(ISOTOPES_GROUP,"/",ISOTOPES_TO_FIND_PREFIX,isamp,"_",1:(length(batches) - 1),sep="")

  #The prefix to WRITE the infos from the hdf5
  feat_input_values <- paste(FEATURES_GROUP,"/",MISSING_FEATURES_VALUES_PREFIX,isamp,"_",1:(length(batches) - 1),sep="")
  isos_input_values <- paste(ISOTOPES_GROUP,"/",ISOTOPES_TO_FIND_VALUES_PREFIX,isamp,"_",1:(length(batches) - 1),sep="")
  feat_input_pos <- paste(FEATURES_GROUP,"/POS_",MISSING_FEATURES_VALUES_PREFIX,isamp,"_",1:(length(batches) - 1),sep="")
  isos_input_pos <- paste(ISOTOPES_GROUP,"/POS_",ISOTOPES_TO_FIND_VALUES_PREFIX,isamp,"_",1:(length(batches) - 1),sep="")
  # if (isamp==24) {
  #   print("")
  # }
  reslist <- extractMissingInformationAllBatches(all_samples[isamp],
                                      all_peaktables[isamp],
                                      isos_input_infos,
                                      feat_input_infos, 
                                      align@rt_correction[[isamp]],
                                      table_iso=isotopes_table, 
                                      dist_c13=dist_iso, 
                                      dm=dm_peaks,
                                      quant=QUANT, 
                                      quant_mean=quant_mean, 
                                      hdf5_file=HDF5_FILE,
                                      extractMissingInformationFun=extractMissingInformation,#For potential parallelization only
                                      margin_mz = margin_dmz,
                                      margin_ppm = margin_ppm,
                                      margin_dmz=margin_dmz,
                                      max_iso = MAX_ISO,
                                      max_charge = MAX_CHARGE,
                                      ppm = PPM,
                                      dmz = DMZ,
                                      detect_isotopes = TRUE)

  #We now write all the possible values in the hdf5 file to avoid fillling th memory in all cases
  for(idx in 1:length(isos_input_infos)){
    infer_idx <- reslist[[2]][[idx]][[1]]
    infer_val <- reslist[[2]][[idx]][[2]]

    #This hsould be written in another way
    isos_idx <- reslist[[1]][[idx]][[1]]
    isos_vals <- reslist[[1]][[idx]][[2]]
    isos_ion_idx <- reslist[[1]][[idx]][[3]]
    if(length(isos_vals)!=0){
      #We write all the isotopes as a single data.frame
      isos_ion_idx <- tryCatch(isos_ion_idx[rep(seq_along(isos_ion_idx),sapply(isos_vals,nrow))],error=function(e){
        print(e)
        print(isos_vals)
        print(isos_idx)
      })
      isos_vals <- do.call(rbind,isos_vals)
      isos_vals$idx <- isos_ion_idx
    }else{
      isos_vals <- data.frame(
        name = EMPTY_ISOTOPES,
        int = 0.0,
        massdiff = 0.0,
        charge = 1,
        idx = 0
      )
    }
    ##WE write the 4 elements to the table
    dummy <- h5write(obj = infer_idx,file = HDF5_FILE,name = feat_input_pos[idx])
    dummy <- h5write(obj = infer_val,file = HDF5_FILE,name = feat_input_values[idx])
    #dummy <- h5write(obj = isos_idx,file = HDF5_FILE,name = isos_input_pos[idx])
    dummy <- h5write(obj = isos_vals,file = HDF5_FILE,name = isos_input_values[idx])
  }
}

# Write the datamatrix by batches in the pekatables

## delete file if exists
output_csv <- gsub('data_prefill_','data_filled_',PATH_DM)
if (file.exists(output_csv)) {
  file.remove(output_csv)
}


###We now perform the actual processing while using the hdf5 peaks to reduce the overhead
for (idx in 1:(length(batches) - 1)) {
  dm <-
    fread(input = PATH_DM,
          nrows = batches[idx + 1] - batches[idx],
          skip = batches[idx],sep="\t")
  colnames(dm) <- cnames

  #Path of the extracted informations into the hdf5
  seq_samples <- 1:num_samples
  feat_input_infos <- paste(FEATURES_GROUP,"/",MISSING_FEATURES_PREFIX,seq_samples,"_",idx,sep="")
  isos_input_infos <- paste(ISOTOPES_GROUP,"/",ISOTOPES_TO_FIND_PREFIX,seq_samples,"_",idx,sep="")
  feat_input_values <- paste(FEATURES_GROUP,"/",MISSING_FEATURES_VALUES_PREFIX,seq_samples,"_",idx,sep="")
  isos_input_values <- paste(ISOTOPES_GROUP,"/",ISOTOPES_TO_FIND_VALUES_PREFIX,seq_samples,"_",idx,sep="")
  feat_input_pos <- paste(FEATURES_GROUP,"/POS_",MISSING_FEATURES_VALUES_PREFIX,seq_samples,"_",idx,sep="")


  # We fill the quantitive informations inferred
  for (iquant in seq_along(quant_cols)) {
    if (is.na(feat_input_infos[iquant])) {
      print(paste('NA found in quant_cols at index',iquant))
      next
    }
    iquant_infer <- h5read(HDF5_FILE,feat_input_infos[iquant])
    if(length(iquant_infer) == 0){
      next
    }
    to_change <- rep(NA,length(iquant_infer))
    sel_col <- quant_cols[iquant]
    infer_pos <- h5read(HDF5_FILE,feat_input_pos[iquant])
    infer_values <- h5read(HDF5_FILE,feat_input_values[iquant])
    to_change[infer_pos] <- infer_values
    dm[[sel_col]][iquant_infer] <- to_change
  }

  # We add the isotopic
  name_col <- rep(NA_character_, nrow(dm))
  dist_col <- rep(NA_real_, nrow(dm))
  abs_col <- rep(NA_real_, nrow(dm))

  for (iquant in 1:length(all_samples)) 
    {
    iquant_isos <- h5read(HDF5_FILE,isos_input_infos[iquant])
    if(length(iquant_isos) == 0)
      next

    isos_values <- h5read(HDF5_FILE,isos_input_values[iquant])

    if(isos_values$name[1]==EMPTY_ISOTOPES){next}

    #We read the indices from the first column of the table
    isos_values <- split(isos_values,isos_values$idx)
    isos_pos <- as.numeric(names(isos_values))

    # names_isos <- rep(NA,length(iquant_isos))
    found_isotopes <- sapply(isos_values, function(x) {
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

    # assemble relative distributions
    found_dists <- sapply(isos_values, function(x) {
      paste(sprintf("%0.4f", x[, "int"] / max(x[, "int"])), collapse = "|")
    })

    #assemble absolute distributions
    found_abs <- sapply(isos_values, function(x) {
      paste(sprintf("%0.0f", x[, "int"]), collapse = "|")
    })

    abs_col[isos_pos] <- found_abs
    name_col[isos_pos] <- found_isotopes
    dist_col[isos_pos] <- found_dists
  }
  
  ###The data matrix is expanded
  infos_idx <- 1:(quant_cols[1] - 1)
  quant_idx <- quant_cols
  new_names <-
    c(colnames(dm)[infos_idx], ISO_NAME_COL, ISO_DIST_COL, ISOABS_DIST_COL, colnames(dm)[quant_idx])
  quant_cols_2 <- quant_cols + 3
  dm <-
    cbind(dm[, ..infos_idx], name_col, dist_col, abs_col, dm[, ..quant_idx])
  colnames(dm) <- new_names
  # we don't write a temp and rename. We just write it to the final destination and keep the prefill around. This avoids a few problems when rerurring the analysis
  if(!file.exists(output_csv)){
    ww <- fwrite(dm, output_csv, sep="\t")
  }else{
    ww <- fwrite(dm, output_csv, sep="\t",append = TRUE)
  }
}
