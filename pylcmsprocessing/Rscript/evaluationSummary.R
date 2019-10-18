suppressWarnings(suppressMessages(library(DBI)))
suppressWarnings(suppressMessages(library(RSQLite)))
suppressWarnings(suppressMessages(library(BiocParallel)))
suppressWarnings(suppressMessages(library(jsonlite)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gghighlight)))
suppressWarnings(suppressMessages(library(stringr)))
suppressWarnings(suppressMessages(library(igraph)))
suppressWarnings(suppressMessages(library(viridis)))
suppressWarnings(suppressMessages(library(ggridges)))
suppressWarnings(suppressMessages(library(lattice)))




####We compute a single metric for each file giving the number of signal detected in num rpelicate-1
###To do se we just ocmpute lal the number of
paths=NULL
computeReproduciblePeak <-  function(eval_table,replicates=NULL,min_sample=3,strict=FALSE){
    ###We read all the name of the evaluation from the peaktable.
  tab <- table(replicates)
  count_replicates <- as.numeric(tab)
  names_replicates <- names(tab)
  ####FOR TESTING PURPOSE
  total_replicates <-  length(names_replicates)

  ###We just the number of rpelicates from the evaluation table
  # tab_replicate <- eval_table[, (ncol(tab) - 2 * num_replicates + 1):(ncol(eval_table) -

  shift <- 0
  if("matched_compound" %in% colnames(eval_table)){
    shift <- 1

  }

  ###Return YES or TRUE depending ofthe function
  to_return <- apply(eval_table,1,function(x,count_replicates,total_replicates,min_sample,shift,strict){
    #####We check if a peak table has been detected



    selsample <- x[(3+shift+1):(3+shift+total_replicates)]
    ####We check if the minumum numbe rof sample is required
    boolpres <-  sapply(seq_along(selsample),function(y,vals,count_replicates,min_thresh){
      if(vals[y]<=min_thresh) return(vals[y]==1)
      return(vals[y]*count_replicates>=(count_replicates-1))

    },vals=selsample,count_replicates=count_replicates,min_thresh=min_sample)
    if(strict){
      return(all(boolpres))
    }else{
      return(any(boolpres))
    }

  },count_replicates=count_replicates,total_replicates=total_replicates,
  min_sample=min_sample,shift=shift,strict=strict)
  return(sum(to_return))
}


####Number of peaks with the detected peak table.
computingStatistic <-  function(eval_table,replicates,min_sample=3){
  ###Number of peak
  if(is.character(eval_table )) eval_table <- read.table(eval_table,sep=",",header=TRUE)
  numpeaks <- nrow(eval_table)
  numreproduciblepeaksallreplicates <- computeReproduciblePeak(eval_table,
                                                            replicates=replicates,min_sample=3,strict=TRUE)
  numreproduciblepeaksreplicates <- computeReproduciblePeak(eval_table,
                                                               replicates=replicates,min_sample=3,strict=FALSE)

  ###Preparatory variable.
  matched <- "matched_compound" %in% colnames(eval_table)
  num_sample <- floor((ncol(eval_table)-2-matched)/3)

  iddetection <- (4+matched+num_sample*0):(3+matched+num_sample*1)
  # idc13 <- (4+matched+num_sample*1):(3+matched+num_sample*2)
  idcv <- (4+matched+num_sample*1):(3+matched+num_sample*2)

  # message("iddetect",iddetection,"idc13",  idc13,"idcv",  idcv)

  mdetection <-  apply(eval_table[,iddetection,drop=FALSE],1,mean,na.rm=TRUE)
  # mc13 <-  apply(eval_table[,idc13,drop=FALSE],1,median,na.rm=TRUE)
  mcv <-  apply(eval_table[,idcv,drop=FALSE],1,mean,na.rm=TRUE)

  ####Untargetted metrics
  npeaks <- nrow(eval_table)


  ###We get the quantile
  intensities <- eval_table[,"mean_intensity"]

  num_quantile <- 5

  qt <- quantile(intensities,probs = seq(0,1,length.out = num_quantile+1))
  qt[1] <- qt[1]-1
  qt[num_quantile+1] <- qt[num_quantile+1]+1
  discrete_lab <- .bincode(intensities,qt)

  ###We calculate the three metrics by quantile
  adetection <- aggregate(mdetection,by=list(discrete_lab),FUN=median,na.rm=TRUE)$x

  # ac13 <- aggregate(mc13,by=list(discrete_lab),FUN=median,na.rm=TRUE)$x
  acv <- aggregate(mcv,by=list(discrete_lab),FUN=median,na.rm=TRUE)$x

  cnames <- c("num_features","reproducible_peaks","always_present","detection_replicate",#"isotope_replicate",
              "cv_replicate",paste("detection_replicate_quantile",1:num_quantile,sep="_"),
              paste("cv_replicate_quantile",1:num_quantile,sep="_"))#paste("isotope_replicate_quantile",1:num_quantile,sep="_"),

  to_return <-  c(total_peaks=numpeaks,
                     reproducible_peaks=numreproduciblepeaksreplicates,
                     always_present=numreproduciblepeaksallreplicates,median(mdetection,na.rm=TRUE),#median(mc13,na.rm=TRUE)
                     median(mcv,na.rm=TRUE),adetection,acv) #ac13

  names(to_return) <- cnames
  return(to_return)
}




###Plot depending og the software on the number of signal versus number of reproducible peak
plotSummary <- function(evaluations, replicates, metrics=NULL,
                        graphical_output = NULL, best =NULL){

  if(is.null(metrics)){
    restable <- sapply(evaluations,computingStatistic,replicates=replicates,simplify = FALSE,USE.NAMES = FALSE)
    metrics <- as.data.frame(do.call(rbind,restable))
  }

  metrics$num_features <- as.numeric(metrics$num_features)
  metrics$reproducible_peaks <- as.numeric(metrics$reproducible_peaks)
  metrics$cv_replicate <- as.numeric(metrics$cv_replicate)


  ####Then we plot the coefficient of variation of the peak versus the number of reporducible peak for each softwre
  ggp <- ggplot(data = metrics,aes(x=reproducible_peaks,y=cv_replicate,color="red"))+
    geom_point(size=2)+xlab("# reproducible peaks")+ylab("# miean CV")+ggtitle("All features : Reproducible Features vs CVs")
  if(!is.null(best)){
    ggp <- ggp +gghighlight(score>=max(metrics$score),use_direct_label = FALSE, use_group_by=FALSE)
  }

  print(ggp)

  ####If the is a targetted column
  if("num_detected_compounds" %in% colnames(metrics)){

    ggp <- ggplot(data = metrics,aes(x=as.numeric(num_reproducible_compounds),y=as.numeric(cv_replicate_compounds),color="red"))+
      geom_point(size=2)+xlab("# reproducible compound peaks")+ylab("# mean compounds CV")+ggtitle("Compounds : Reproducible detected Features vs CVs")
    if(!is.null(best)){
      ggp <- ggp +gghighlight(score>=max(metrics$score),use_direct_label = FALSE, use_group_by=FALSE)
    }

    print(ggp)

    ggp <- ggplot(data = metrics,aes(x=as.numeric(num_detected_compounds),y=as.numeric(cv_replicate_compounds),color=software,shape=software))+
      geom_point(size=2)+xlab("# compound peaks")+ylab("# mean compounds CV")+ggtitle("Compounds : Detected Features vs CVs")
    if(!is.null(best)){
      ggp <- ggp +gghighlight(score>=max(metrics$score),use_direct_label = FALSE, use_group_by=FALSE)
    }

    print(ggp)
  }
  ###We close the graphical device.

  return(metrics)
}

score_func_targetted <-  function(x,fix=0.5){
  as.numeric(x["num_reproducible_compounds"])*(fix-as.numeric(x["cv_replicate_compounds"]))
}

score_func_untargetted <-  function(x,fix=0.5){
  as.numeric(x["reproducible_peaks"])*(fix-as.numeric(x["cv_replicate"]))
}

####Returnt he index of the best perming peak picking depending of the mass spectrometer

findBestPeakPicking <- function(metrics,type=c("untargetted","targetted"),score=NULL){
  type <-  match.arg(type)
  if(is.null(score)){
    if(type=="targetted"){
      score <- score_func_targetted
    }else if(type=="untargetted"){
      score <- score_func_untargetted
    }
  }else if(typeof("score")!= "closure"){
    stop("score should be a closure.")
  }
  vscore <- apply(metrics,1,score)
  pmax <- which.max(vscore)
  list(pmax,vscore)
}



plotCVIntensity <-  function(eval_table,vtitle="Best parameter sets",num_quantile=10){
  pcv <- str_detect(colnames(eval_table),fixed("CV_replicate"))
  fullcv <- apply(eval_table[,pcv,drop=FALSE],1,mean,na.rm=TRUE)
  psel <- !sapply(fullcv,is.nan)

  lint <- log10(eval_table[psel,"median_intensity"])

  qt <- quantile(lint,probs = seq(0,1,length.out = num_quantile+1))
  qt[1] <- qt[1]-1
  qt[num_quantile] <- qt[num_quantile]+1
  discrete_lab <- .bincode(lint,qt)

  df <- data.frame(cv=fullcv[psel],logint=lint, intq=discrete_lab)

  intlab <- paste(sprintf("%.2e",qt[1:(length(qt)-1)]),"-\n",sprintf("%.2e",qt[2:length(qt)]),sep="")

  ggp <- ggplot()+  geom_boxplot(data=df,aes(x=as.factor(intq),y=cv,fill=as.factor(intq)))+xlab("Intensity quantile")+
    ylab("median CV")+ggtitle(vtitle)+
    scale_fill_manual(name = "Intensity Q", labels = intlab,values = viridis(length(intlab)))
  print(ggp)
}

# heatMapCompounds <- function(datamatrix,matched,id,samples=NULL,compounds=NULL){
#
#   p_matched <- which(!is.na(matched))
#
#   # matched <-  which(!is.na(matched))
#   if(length(matched)==0) stop("impossible to match anything")
#
#   datamatrix
#
#   ###WE tak ethe log of the dataset
#   pna <- which(is.na(datamatrix),arr.ind = TRUE)
#
#   datamatrix[pna] <- 1
#
#   transpf <- as.matrix(log10(datamatrix[3:ncol(datamatrix)]))
#
#   found_compounds <- !is.na(matched)
#   pfound <-  matched[p_matched]
#
#   if(is.null(samples)){
#     samples <-  paste("sample",1:ncol(transpf),sep=" ")
#   }
#   if(is.null(compounds)){
#     compounds <-  paste("compounds",1:length(pfound),sep=" ")
#   }else{
#     compounds <- compounds[p_matched]
#   }
#   ###We add this thing of the data
#   data <- expand.grid(X=samples, Y=compounds)
#   data$Z <- as.numeric(transpf[pfound,,drop=FALSE])
#   levelplot(Z ~ X*Y, data=data  , xlab="Samples" ,scales=list(x=list(cex=.3),y=list(cex=.7)),
#             ylab="Compounds", col.regions = viridis_pal()(100)   , main="")
#   message("levelplot done")
# }



evaluatepeakpicking <-  function(peakpickings,replicates,best_param_path,path_graphics=NULL,vmetrics=NULL){

  if(!is.null(path_graphics)){
    pdf(path_graphics)
  }

   evaluations <- peakpickings$evaluation

  if(is.null(vmetrics)){
    vmetrics <- plotSummary(evaluations,replicates)
  }

  lbpf <- findBestPeakPicking(vmetrics,type="untargetted")
  bpf <- lbpf[[1]]

  vmetrics$score <- lbpf[[2]]

  peaktables <- peakpickings$peaktable
  bestpf <- peaktables[bpf[1]]

  message("best peak picking is ",bestpf)

  best_peaktable <- read.table(bestpf,sep=",",header=TRUE)

  best_eval <- read.table(evaluations[bpf[1]],sep=",",header=TRUE)

  plotSummary(evaluations,replicates)

  plotCVIntensity(best_eval,vtitle=basename(peakpickings$parameter[bpf]),num_quantile=10)

  if(!is.null(path_graphics)){
    dev.off()
  }

  ####We copy the best parameters files.
  file.copy(peakpickings$parameter[bpf[1]],best_param_path,overwrite=TRUE)

  return(invisible(list(metrics=vmetrics,best=bpf[1])))
}

args <- commandArgs(trailingOnly = TRUE)

PATHDB <- args[1]
FIGURE <- args[2]
BEST_PARAM <- args[3]
OUTPEAKTABLES <- args[4]
db <- dbConnect(RSQLite::SQLite(),PATHDB)
ppts <- dbGetQuery(db,"SELECT * FROM peakpicking")
replicates <- dbGetQuery(db,"SELECT replicate FROM samples")[,1]
dbDisconnect(db)

res <- evaluatepeakpicking(ppts,replicates,best_param_path=BEST_PARAM,path_graphics=FIGURE)
file_res_n <- file.path(dirname(args[1]),"metrics_eval.csv")
write.table(res$metrics,file=file_res_n,row.names=FALSE,sep=",")


db <- dbConnect(RSQLite::SQLite(),PATHDB)
ppts2 <- dbGetQuery(db,paste("SELECT output FROM processing WHERE peakpicking=",ppts$id[res$best],sep=""))[,1]
peaktable <- dbGetQuery(db,paste("SELECT peaktable FROM peakpicking WHERE id=",ppts$id[res$best],sep=""))[,1]
dbDisconnect(db)
####We copy all of these in the new directory
for(p in ppts2){
  file.copy(p,file.path(OUTPEAKTABLES,basename(p)),overwrite=TRUE)
}
file.copy(peaktable,file.path(OUTPEAKTABLES,basename(peaktable)),overwrite=TRUE)
