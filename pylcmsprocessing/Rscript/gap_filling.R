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
MAX_FILES <- 10

sink(file=stdout())

##Argument passed by Python
args <- commandArgs(trailingOnly = TRUE)

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
NUM_FILES <- as.integer(args[12])
NUM_WORKERS <- as.integer(args[13])

if(NUM_FILES>MAX_FILES) NUM_FILES <- MAX_FILES

TEMP_NAME <- file.path(dirname(PATH_FILLED),"temp_transfer.csv")

###This is jsut for evaluation
TEMP_FILLED <- paste(TEMP_NAME, "non_filled.csv", sep = "_")
nullvar <- file.copy(PATH_DM, TEMP_FILLED)
print(TEMP_NAME)
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

all_infos <-
  dbGetQuery(
    dbb,
    "SELECT path,output_ms,types FROM samples INNER JOIN processing on samples.id=processing.sample WHERE level='MS1' AND output_ms!='NOT PROCESSED' AND valid=1"
  )

all_peaktables <- all_infos[,2]
all_samples <- all_infos[,1]
# if(FALSE){
#   all_peaktables <- file.path("U:/processing/out/sammy slaw 2/CENTWAVE/peaktables",basename(all_peaktables))
#   all_samples <- file.path("U:/processing/out/sammy slaw 2/mzMLs-delete-when-done",basename(all_samples))
# }

optim_idx <- seq_along(all_samples)

if("QC" %in% all_infos[,3]){
  optim_idx <- which(all_infos[,3]=="QC")
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
dm <- fread(PATH_DM, sep = "\t", select = c("mz", "rt"))
max_feat <- nrow(dm)
BY_BATCH <- 2000

rm(dm)
batches <- seq(1, max_feat, by = BY_BATCH)
batches[1] <- 0
if (batches[length(batches)] != max_feat) {
  batches <- c(batches, max_feat + 2)
} else{
  batches[length(batches)] <- max_feat + 2
}


###We define the gap filling and isotopes detection fucntion

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
           mean_int,
           margin_dmz = 0.003,
           margin_ppm = 10,
           margin_mz = 0.001,
           max_iso = 4,
           max_charge = 2,
           ppm = 8,
           dmz = 0.002,
           detect_isotopes = TRUE){
    ###We create default vector to return if needed
    def_isotopes <- rep(NA, length(isotopes))
    def_missing <- rep(0.0, length(to_infer))
    suppressWarnings(suppressMessages(library(pracma, warn.conflicts = FALSE)))
    suppressWarnings(suppressMessages(library(data.table, warn.conflicts = FALSE)))
    
    xraw <- praw
    if(is.character(praw)){
      suppressWarnings(suppressMessages(library(xcms, warn.conflicts = FALSE)))
      xraw <- suppressWarnings(suppressMessages(xcmsRaw(praw)))
    }
    peaks <- fread(peaks)
    
    noise_level <- quantile(peaks$intensity, probs = 0.03)
    ort <- order(dm[infer, ][[4]])
    ###small cheat the correction is the closest
    all_corrs <- rep(0,length(infer[ort]))
    if(is.na(align[1, 1])){
      dev_rt <- rep(0,length(infer[ort]))
    }else{
      vmatch  <-
        xcms:::findEqualGreaterM(align[, 1], values = dm[infer[ort], ][[4]])
      align <- rbind(align, align[nrow(align), ])
      dev_rt <- align[vmatch, 2]
    }
    ###integrate the missing peak
    rt_corr <- dm[[4]][infer[ort]] + dev_rt
    
    ###we first perform gap filling
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
        
        integ_fun <- function(x,y){
          trapz(x,y)
        }
        if (quant != "intensity") {
          integ_fun <- function(x,y){
            max(y)
          }
        }
        
        rt_min <- rt - peakwidth / 2
        rt_max <- rt + peakwidth / 2
        #We shift the peak if it is negativer
        if(rt_min<0){
          rt_max <- rt_max+abs(rt_min)
          rt_min <- 0
        }
        
        
        mean_mz <- (min_mz+max_mz)/2
        
        ###Margin in ppm.
        margin_mz <- max(margin_dmz,mean_mz*margin_ppm*1e-6)
        min_mz <- min_mz - margin_mz
        max_mz <- max_mz + margin_mz
        total_extension <- margin_mz
        # cat("mzr:",paste(sprintf("%0.4f",c(min_mz,max_mz)),collapse = "-"))
        
        tval <-
          tryCatch(
            rawEIC(
              xraw,
              rtrange = c(rt_min * 60, rt_max * 60),
              mzrange = c(min_mz, max_mz)
            ),
            error = function(e) {
              # print(c(rt_min * 60, rt_max * 60))
              # print(c(min_mz, max_mz))
              # print(e)
              return(NA)
            }
          )
        if (length(tval) == 1 && is.na(tval))
          return(0)
        val <- 0
        ##We extend the limit one more time.
        max_extension <- 0.07
        current_extension <- 0
        old_int <- -1
        current_int <- integ_fun(x = xraw@scantime[tval[[1]]], tval[[2]])
        # if(abs(min_mz-1397.37)<0.05 & abs(rt-1.41)<0.15) browser()
        ###We keep expanding the trace mass.
        while ((abs(expected_intensity-current_int)<=abs(expected_intensity-old_int))&
               (old_int != current_int) & (current_int<expected_intensity)&
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
          total_extension <- total_extension+margin_mz
        }
        val <- old_int
        if(val==-1) val<- 0
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
          mean_int[infer[ort]],
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
        }
      )
    
    inferred_values <-
      inferred_values[order(ort, decreasing = FALSE)]
    
    visotopes <- NULL
    if(detect_isotopes){
      
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
    }
    return(list(inferred_values, visotopes))
  }




# Before doing any batches we try to optimize filling the mz_filling windows as a function
# of the parameters of everything
dm <-
  fread(input = PATH_DM,sep="\t")
colnames(dm) <- cnames

sel_columns <-
  c("min_mz", "max_mz", "mean_mz", "mean_rt", "mean_peakwidth")
if (max(dm[["mean_peakwidth"]]) > max(dm[["mean_rt"]])) {
  dm[["mean_peakwidth"]] <- dm[["mean_peakwidth"]] / 60
}
dm_peaks <- dm[, ..sel_columns]

quant_prefix <-
  paste(str_split(cnames[length(cnames)], fixed("_"))[[1]][1], "_", sep =
          "")

quant_cols <- which(startsWith(cnames, quant_prefix))

align <- readRDS(PATH_MODEL)
quant_cols_optim <- quant_cols[optim_idx]

isotopes_to_extract <-
  vector(mode = "list", length = length(quant_cols))
for (ic in seq_along(isotopes_to_extract)) {
  isotopes_to_extract[[ic]] <- numeric(0)
}
dist_iso <- NULL
# message("Inferring",all_samples,"vs",,)
###We also extract the data to infer
to_infer <-
  apply(dm[, ..quant_cols_optim], 2, function(x) {
    which(is.na(x))
  })


optim_intensity <- apply(dm[,..quant_cols_optim],1,mean,na.rm=TRUE)



###We keep increasing the tolerance until the distribution stat to match.
optimizeParameters <-function(praws,peaks,isotopes,infer,
                             align,table_iso,dist_c13,
                             dm,quant,optim_intensity,
                             max_iso = 4,max_charge = 2,
                             ppm = 8,dmz = 0.002, by_file = 300){
  
  ###We select at most 100 values by files.
  infer <- lapply(infer,FUN = function(x,max_size){
    if(length(x)>max_size){
      return(sample(x,max_size))
    }
    return(x)
  },max_size=by_file)
  
  
  suppressWarnings(suppressMessages(library(xcms, warn.conflicts = FALSE)))
  ###Reading the file only once
  xraws <- suppressWarnings(suppressMessages(sapply(praws,xcmsRaw)))

  ###One ce the best margin has been found we can recompute the distribution
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
                     infers,aligns, FUN = extractMissingInformations,
                     MoreArgs=list(table_iso=table_iso,dist_c13=dist_c13,dm=dm,quant=quant,mean_int=optim_intensity,
                                   margin_mz = x,
                                   margin_ppm = x[2],
                                   margin_dmz = x[1],
                                   max_iso = max_iso,max_charge = max_charge,ppm = ppm,
                                   dmz = dmz,detect_isotopes = FALSE),SIMPLIFY = FALSE)
    
    score <- 0
    for(idx in seq_along(xraws)){
      to_average <- abs(log10(vcomps[[idx]][[1]]+2)-log10(optim_intensity[infers[[idx]]]+2))
      score <- score + mean(to_average[optim_intensity[infers[[idx]]]!=0],na.rm = TRUE)
    }
    score <- score/length(xraws)
    return(score)
  }


  #This is potentially very costly in memory so we reduce it
  dmz_seq <- exp(seq(log(0.0005),log(0.05),length=5))
  ppm_seq <- exp(seq(log(2),log(30),length=5))
  grid_seq <- expand.grid(dmz_seq,ppm_seq)
  veval <- apply(grid_seq,1,FUN = score_gap_filling,xraws=xraws,
        peaks=peaks,isotopes=isotopes,infer=infer,
        align=align,table_iso=table_iso,dist_c13=dist_c13,
        dm=dm,quant=quant,optim_intensity=optim_intensity,max_iso = max_iso,
        max_charge = max_charge,ppm = ppm,
        dmz = dmz)
  
  best_par <- unlist(grid_seq[which.min(veval),])
  best_val <- min(veval)


  ###one we found a local minima, for the three minimal value we just od a gradient descent.
  dmz_step <- 0.0001
  ppm_step <- 0.2
  new_val <- -1
  new_par <- unlist(best_par)
  MAX_ITS <- 10
  current_it <- 1
  # print(best_par)
  # print(best_val)
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
    
    # cat("val_dmz:",val_dmz)
    # print("best_val:",best_val)
    # print("dmz_step:",dmz_step)
    
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
    
    # cat(currenti_it,new_par,new_val)
    current_it <- current_it+1
  }
  return(c(best_par,best_val))
}

# cat("Optimizing parameters")

##We only keep the best parameters
temp_par <- tryCatch(optimizeParameters(all_samples[optim_idx],all_peaktables[optim_idx],isotopes_to_extract,
                               to_infer,align@rt_correction,
                               table_iso= isotopes_table,dist_c13 = dist_iso,
                               dm = dm_peaks,quant = QUANT,optim_intensity = optim_intensity,
                               max_iso = 4,max_charge = 2,
                               ppm = 8,dmz = 0.002),error=function(e){return(c(0.005,10))})

# cat("Optimization finished.")
margin_dmz <- temp_par[1]
margin_ppm <- temp_par[2]
rm(dm)



###We now perfomr the actual processing
for (idx in 1:(length(batches) - 1)) {
  dm <-
    fread(input = PATH_DM,
          nrows = batches[idx + 1] - batches[idx] - 1,
          skip = batches[idx],sep="\t")
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
  
  ###We also extract the data to infer
  to_infer <-
    apply(dm[, ..quant_cols], 2, function(x) {
      which(is.na(x))
    })
  

  expected_intensity <- apply(dm[,..quant_cols],1,mean,na.rm=TRUE)
  
  ###We first determine the gap-filling parameters on a set of QCs file.
  # cat("margin_ppm",sprintf("%0.2f",margin_ppm),"\n")
  # cat("margin_dmz",sprintf("%0.5f",margin_dmz),"\n")
  # 
  # for(idx in seq_along(all_samples)){
  # 
  #   exm <- extractMissingInformations(   all_samples[idx],
  #                                 all_peaktables[idx],
  #                                 isotopes_to_extract[[idx]],
  #                                 to_infer[[idx]],
  #                                 align@rt_correction[[idx]],
  #                                 dm = dm_peaks,
  #                                 quant = QUANT,
  #                                 table_iso = isotopes_table,
  #                                 dist_c13 =
  #                                   dist_iso,
  #                                 margin_mz = margin_dmz,
  #                                 margin_ppm = margin_ppm,
  #                                 margin_dmz=margin_dmz,
  #                                 max_iso = MAX_ISO,
  #                                 max_charge = MAX_CHARGE,
  #                                 mean_int = expected_intensity,
  #                                 ppm = PPM,
  #                                 dmz = DMZ)
  # 
  # }


  
  vmap <-
    bpmapply(
      FUN = extractMissingInformations,
      all_samples,
      all_peaktables,
      isotopes_to_extract,
      to_infer,
      align@rt_correction,
      MoreArgs = list(
        dm = dm_peaks,
        quant = QUANT,
        table_iso = isotopes_table,
        dist_c13 =
          dist_iso,
        margin_mz = margin_dmz,
        margin_ppm = margin_ppm,
        margin_dmz=margin_dmz,
        max_iso = MAX_ISO,
        max_charge = MAX_CHARGE,
        mean_int = expected_intensity,
        ppm = PPM,
        dmz = DMZ
      ),
      SIMPLIFY = FALSE,
      BPPARAM = bpp
    )
  # }
  ###We fill the column
  for (iquant in seq_along(quant_cols)) {
    if (length(to_infer[[iquant]]) == 0){
      cat("No imputed values.")
      next
    }
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
  if(!file.exists(TEMP_NAME)){
    ww <- fwrite(dm, TEMP_NAME,sep="\t")
  }else{
    ww <- fwrite(dm, TEMP_NAME,sep="\t",append = TRUE)
  }
}
ww <- file.rename(PATH_DM, TEMP_FILLED)
ww <- file.rename(TEMP_NAME, PATH_DM)

# ww <- file.rename(TEMP_NAME, PATH_DM)
