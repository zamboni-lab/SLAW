####
suppressMessages(library(DBI))
suppressMessages(library(RSQLite))
suppressMessages(library(igraph))
suppressMessages(library(maxmatching))
suppressMessages(library(proFIA))
suppressMessages(library(stringr))
suppressMessages(library(BiocParallel))
suppressMessages(library(ggplot2))
suppressMessages(library(gghighlight))
suppressMessages(library(viridis))


####Compute an evaluation summary for a given data table

###For each peaktable, calculate a single of each of the performance metrics
###for each set of replicates for each variable.
getEvaluationMatrix <- function(peaktable, output, replicates) {

  suppressMessages(library(DBI))
  suppressMessages(library(RSQLite))
  suppressMessages(library(igraph))
  suppressMessages(library(maxmatching))
  suppressMessages(library(proFIA))
  suppressMessages(library(ggplot2))
  suppressMessages(library(gghighlight))
  suppressMessages(library(stringr))


  evaluateDetection <-
    function(vmatrix, replicates) {

      freqt <-  table(replicates)
      matDetection <- matrix(0.0, nrow = nrow(vmatrix), ncol = length(freqt))

      ###We collect the precentage
      allAggregate <-  apply(vmatrix, 1, function(row, vreplicate) {
        aggregate(
          x = row,
          by = list(vreplicate),
          FUN = function(x) {
            temp <- sum(!is.na(x))
            return(c(temp, temp / length(x)))
          }
        )
      }, vreplicate = replicates)

      cnames <- paste("Detection_frequency_replicate_",names(freqt),sep="")
      ####We return this number of rall the signal
      vres <- do.call(rbind,sapply(allAggregate, function(x) {

        return((x$x[,2:ncol(x$x),drop=FALSE]))
      },simplify = FALSE))
      if(ncol(vres)>nrow(vres))
      vres <-  t(vres)
      colnames(vres) <- cnames

      ####We add columns names
      return(vres)
    }


    ###Evaluation of detection reproducibility
  evaluateFractionDetectionReplicate <-  function( peaktable, replicates) {
    if (is.character(peaktable)) {
      peaktable <- read.table(peaktable, sep = ",", header = TRUE)
    }
    return(evaluateDetection(peaktable, replicates))
  }



  evaluateCV <-
    function(vmatrix, replicates) {
      freqt <-  table(replicates)
      ###We collect the precentage
      allAggregate <-  apply(vmatrix, 1, function(row, vreplicate) {
        aggregate(
          x = row,
          by = list(vreplicate),
          FUN = function(x) {
            ###We calculate the CV
            pok <- which(!is.na(x))
            if(length(pok)>=2){
              return(sd(x[pok])/mean(x[pok]))
            }else{
              return(NA)
            }
          }
        )
      }, vreplicate = replicates)
      cnames <- paste("CV_replicate_",names(freqt),sep="")
      ####We return this number of rall the signal
      vres <- do.call(rbind,sapply(allAggregate, function(x) {
        return((x$x))
      },simplify=FALSE))
      if(ncol(vres)>nrow(vres))
      vres <-  t(vres)
      colnames(vres) <- cnames
      ####We add columns names
      return(vres)
    }

  evaluateCVReplicate <-  function(peaktable, replicates) {
    if (is.character(peaktable)) {
      peaktable <- read.table(peaktable, sep = ",", header = TRUE)
    }
    return(evaluateCV(peaktable, replicates))
  }


  ###GIven a data matrix  Compute the fraction of detected isotope peak
  findC13Replicate <-  function(vmatrix,replicates,
                                ppm = 15,
                                dmz = 0.003,
                                tolrt = 20,
                                frac = 0.8) {
    ###We find the adapted documentary.
    C13mzdiff <- 1.003355
    mzC13 <-  vmatrix[, "mz"] + C13mzdiff



    urepli <-  unique(replicates)


    matIso <-  matrix(NA_real_,nrow=nrow(vmatrix),ncol=length(urepli))

    ###We get the position of the possible isotopic peaks
    vm <-  proFIA:::fastMatchPpm(mzC13, vmatrix[, "mz"], ppm = ppm)#, dmz = dmz)
    posnf <- sapply(vm, is.null)
    if(sum(posnf)!=0){
      for(pi in which(posnf)){
        vm[[pi]] <- numeric(0)
      }
    }


    posf <-  !posnf
    ####We filter them by retnetion time now
    vm[posf] <- mapply(vm[posf], which(posf), FUN=function(x, xi, rts, tolrt) {
      x[abs(rts[xi] - rts[x]) < tolrt]
    }, MoreArgs = list(rts = vmatrix[, "rt"], tolrt = tolrt),SIMPLIFY = FALSE)


    ###We evalaute the intensities ratio eventually
    # sample(seq(1, ncol(vmatrix)), size = ncol(vmatrix))

    posnf <- sapply(vm, length) != 0
    # mCH2 <-  14.01565
    # mCH3 <-  15.02348
    # pC13 <- 1.109e-2


    ####We determine the noise level in each sample 3/100 of detected peaks.
    vnoise <- apply(vmatrix[, 3:ncol(vmatrix),drop=FALSE], 2,FUN =  function(x) {
      quantile(x[x != 0], 3 / 100,na.rm = TRUE)
    })
    ###we can calculate for each compound its observability given the level of noise in the sample

    ###NA_real_ mean that the data are too noisy an dhtis peak in unspottable anyway.
    ####We eventually filter by intensity depending of the isotopoes
    # message(paste(dim(vmatrix),collapse = "|"))
    tres <-
      mapply(vm[posf], which(posf), FUN = function(x, xi, mz, vmat, vnoise, repli) {

        ####We first caculate if the isotopic peak si supposed ot be observable
        mass <- mz[xi]
        maxC <- floor((mz[xi] - 2 * 15.02348) / 14.01565) + 2
        ####We calculate the theoric highest intensity for sample
        # message(paste(dim(vmat),collapse = "|"))
        theoric_highest_intensity <- maxC * 1.109e-2 * vmat[xi, 3:ncol(vmat)]
        ###Now we split the results on rpelicates
        ###Only the firs tline is considered in every case
        # message("Comparing Mono : mz ",mz[xi]," C13 : mz ; ",mz[x[1]],"\n")
        gagg <- aggregate(1:nrow(repli),by=list(repli$replicate),FUN=function(xidx,vint_mono,vint_C13,vnoise,maxC){

          ####We keep a multiplicative factor corresponding to the number of na.
          pok <- which(!is.na(vint_mono[xidx]))
          if(length(pok)==0){
            return(NA_real_)
          }

          ###some na will be induced if the peak is not detected.
          theoric_highest_intensity <- maxC * 1.109e-2 * vint_mono[xidx[pok]]

          ###Then the process is done for each variable and replicate.
          pnoise <- (theoric_highest_intensity < vnoise[xidx[pok]])

          ###Impossible ot spot
          if (all(pnoise)) {
            return(NA_real_)
          }


          ###We check which of the peaks has a coherent intensity
          possible <-  !pnoise

          ####They are possible only if they are smaller than the theoric max intnesity
          correct <-  possible & (vint_C13[xidx[pok]] < theoric_highest_intensity)
          # message(paste(correct,sep="|"))

          # message("Comparing Mono : int ",paste(sprintf("%0.1f",vint_mono[xidx]),collapse = "|"),
          #         " C13 : int : ",paste(sprintf("%0.1f",vint_C13[xidx]),collapse = "|"), " theoint ",
          #         paste(sprintf("%0.1f",theoric_highest_intensity),collapse = "|")," final ",sum(correct,na.rm=TRUE)/length(pok),"\n")
          #

          ###Now we return the fraction of the possible on the total length
          return(sum(correct,na.rm = TRUE)/length(pok))

        },vint_mono=vmat[xi, 3:ncol(vmat)],vint_C13=vmat[x[1],3:ncol(vmat)],vnoise=vnoise,maxC=maxC)
        ###We return a vector containing all the values.
        gagg$x
      },MoreArgs = list(mz=vmatrix[,"mz"], vmat=vmatrix, vnoise=vnoise, repli=replicates, frac = 0.8),SIMPLIFY = FALSE)
    cnames <- paste("Isotopes_frequency_replicates_",urepli,sep="")
    matIso[posf,] <- do.call("rbind",tres)
    colnames(matIso) <- cnames
    return(matIso)
  }

  ###we get the data matrix
  detection <- evaluateFractionDetectionReplicate(peaktable[,3:ncol(peaktable)],replicates)
  CV <- evaluateCVReplicate(peaktable[,3:ncol(peaktable)],replicates)
  # valC13 <- findC13Replicate(peaktable[,3:ncol(peaktable)],replicates)
  ###We get an estination of th eintensity, for each components
  vintensity <- apply(peaktable[,1:ncol(peaktable)],1,mean,na.rm=TRUE)

  df_intensity <-  data.frame(mean_intensity=vintensity)


  ddmat <- cbind(peaktable[,c("mz","rt")],df_intensity,detection,CV)

  ####Now we extract supplemetary for the carry

  write.table(x = ddmat,file = output,sep=',',row.names = FALSE)
}



args <- commandArgs(trailingOnly = TRUE)
DATAMATRIX <- args[1]
REPLICATES <- args[2]
OUTPUT <- args[3]


message("evaluating  data matrix ",DATAMATRIX)


peaktable <- read.table(DATAMATRIX,sep=",",header=TRUE)
fopen <- file(REPLICATES,"r")
replicates <- readLines(fopen)
close.connection(fopen)

getEvaluationMatrix(peaktable, OUTPUT, replicates)
