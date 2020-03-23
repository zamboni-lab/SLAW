#' @include classes.R
#' @include correctRTs.R
#' @include utils.R

# setClass("LCMSAlignerModel",slots=list(references="LCMSAlignerReferences",
#                                        peaks="data.frame",
#                                        data="list",
#                                        files="character",
#                                        num_files="numeric",
#                                        rt_correction="list"))


#' Initialisation ofhte LCMSAlignerModel
#'
#' @param peaktables A list of data.frame containing everything
#' @param ref_file The reference file to use a the first stpe fo the dat matrix
#' @param num_peaks Number of reference peaks to consider
#' @param col_int The name of the column with intensity
#' @param ... Supplemetary arguements which will be passed to the nez dataset.
#'
#' @return
#' @export
#'
#' @examples
LCMSAlignerModel <-
  function(paths,
           output,
           ref_file = NULL,
           num_peaks = 150,
           save_interval = 5,
           col_int = "intensity",
           lar = NULL,
           rt = 1,
           rt_dens = 0.05,
           ppm = 12,
           dmz = 0.007,
           int = 1,
           sep_table = ",",
           alignment_threshold = 8,
           algorithm=c("density","nn"),
           graphical = FALSE,
           bpp=NULL,
           seed = 512,supp_data=character(0),
           ...) {
    if (!is.character(paths)) {
      stop("Paths should be a character vector containi gthe names of the pekabtales")
    }

    if(is.null(bpp)){
      bpp <- bpparam()
    }

    algorithm <- match.arg(algorithm)

    if ((length(paths) == 1) && file.exists(paths)) {
      ###Then it s anobject whihc need ot bextracted.
      lam <- readAlignment(paths)
      return(lam)
    }

    peaktables <-
      sapply(
        paths,
        read.table,
        sep = sep_table,
        header = TRUE,
        simplify = FALSE
      )
    ###We perform some security check
    if (length(peaktables) > 20) {
      stop("No more than 20 peaktables as references for computational purpose.")
    }

    lam <- new("LCMSAlignerModel")
    if (!is.null(lar)) {
      message("lar furnished num_peaks argument is ignored.")
      lam@references <- lar
    } else{
      lam@references <-
        LCMSAlignerReferences(
          peaktables,
          num_peaks = num_peaks,
          col_int = col_int,
          seed=seed,
          graphical = graphical,
          bpp=bpp,
          ...
        )
    }

    lam@storage <-  LCMSAlignerDataStorage(output)
    lam@save_interval <- save_interval

    ####We tune the parameters eventually
    params(lam@references) <-
      list(
        rt = rt,
        rt_dens = rt_dens,
        ppm = ppm,
        dmz = dmz,
        int = int,
        sep_table = sep_table,
        alignment_threshold = alignment_threshold,
        aligner_type = algorithm
      )

    lam@max_id <- 0

    ###Parameter settled
    p_ref <- NULL
    if (is.null(ref_file)) {
      p_ref <- which.max(sapply(peaktables, nrow))
      ref_file <- peaktables[[which.max(sapply(peaktables, nrow))]]
    } else{
      p_ref <- match(ref_file, paths)
      if (is.na(p_ref))
        stop('ref_file should be a name corresponding to the reference file used.')
    }

    ref_file <-
      fortifyPeaktable(peaktables[[p_ref]], lar = lam@references)

    ####WWe generate the new object

    lam@files <- LCMSAlignerFileHandler()

    set.seed(seed)



    ###Alll the element of the peak tbale are designated as a dataset
    lam@peaks <-
      data.frame(
        mz = numeric(0),
        rt = numeric(0),
        rtdev = numeric(0),
        rtfactor = numeric(0),
        num = numeric(0),
        current_num = numeric(0),
        id = numeric(0),
        stringsAsFactors = FALSE
      )
    lam@data <- matrix(ncol = 6+length(supp_data), nrow = 0)
    colnames(lam@data) <- c('mz', 'rt', 'rt_cor', 'int',supp_data, 'sample', 'id')

    # lam <- alignPeaktable(lam,
    #                       paths[p_ref])

    return(lam)
  }


#' Align all the file form a directory to an exisiting LCMSAlignerModel eventually
#'
#' @param directory The directory to monitor
#' @param path_model The path to store the model
#' @param output The directory of output for for the ocnstructed block file
#' @param save_interval The number of aligned peak tbale between save
#' @param num_file If A new modle is created, how many files are used to find the references peaks
#' @param reset Shall the processing be restarted from the beginning eventually
#' @param span Smoothing parameters contorlling the loess
#' @param n_clusters Numbers of clusters evnetually tested.
#' @param ... Arguments which will be passed to LCMSAlignerModel
#'
#' @return The LCMSAlignerObjects
#' @export
#'
#' @examples
#' print("Examples to be put here")
LCMSAlignerModelFromDirectory <-
  function(directory,
           path_model,
           output,
           save_interval = 5,
           num_file = 10,
           reset = FALSE,
           threshold = 5,
           span = 0.6,
           n_clusters = 15,
           supp_data=c("peakwidth","SN","right_on_left_assymetry"),
           ransac_l1=0.05,
           ransac_dist_threshold=0.1,
           maxAlign = 100,
           seed=512,
           max_cor = 0.1,
           clustering=TRUE,
           bpp=NULL,
           graphical=FALSE,
           ...) {
    if(length(directory)==1){
    message("Processing of directory ", directory, ".")
    }else{
      message("Processing of ", length(directory), " files.")
    }
    if(is.null(maxAlign)) maxAlign <- 100000
    if(is.null(bpp)) bpp <- bpparam()

    ###Two cases or the furnished vector is a set of chracter or the furnished data table is a set of character
    online <- TRUE
    if((length(directory)==1) && dir.exists(directory)){
        all_files <- list.files(directory, full.names = TRUE)
    }else{ ###It is a set fo files eventually.
        all_files <- directory
        online <- FALSE
    }

    lam <- NULL
    ###If the file is null we intialize it.
    if (file.exists(path_model)) {
      lam <-
        LCMSAlignerModel(
          paths = path_model,
          save_interval = save_interval,
          output = output,
          n_clusters = n_clusters,
          supp_data=supp_data,
          graphical = graphical,
          bpp = bpp,
          ...
        )
      if (reset) {
        lam <- restart(lam,supp_data=supp_data)
      }
    } else{

      ###Change oa pradigm we tak the biggest ifle
      vsize <- file.info(all_files)$size
      ref_files <- all_files[order(vsize,decreasing=TRUE)[1:min(num_file, length(all_files))]]

      lam <-
        LCMSAlignerModel(
          paths = ref_files,
          save_interval = save_interval,
          output = output,
          n_clusters = n_clusters,
          supp_data=supp_data,
          graphical = graphical,
          seed=seed,
          bpp = bpp,
          ...
        )
    }
    ###We check the file to be processed
    vpro <- isProcessed(lam, all_files)
    all_files <- all_files[!vpro]
    num_align <- 0

    while (length(all_files) != 0) {
      num_align <- num_align+1
      if(num_align>maxAlign) break
      if(online){
        all_files <- list.files(directory, full.names = TRUE)
      }
      vpro <- isProcessed(lam, all_files)
      all_files <- all_files[!vpro]
      if (length(all_files) != 0) {
        message("Found ", length(all_files), " files to process.")
      }
      for (nf in all_files) {
        ###We align each file individually
        lam <- alignPeaktable(lam,
                              nf, path_aligner = path_model, span = span,
                              supp_data=supp_data,ransac_l1=ransac_l1,
                              ransac_dist_threshold = ransac_dist_threshold,
                              lim = max_cor, graphical = graphical, bpp = bpp)
      }
      ## we initilaize the clusteritng
      lam@clustering <- 1:nrow(lam@peaks)
    }
    ###We initialize the

    ##We do the ifnal clustering mode
    if(nrow(lam@data)!=0){
      lam <- saveAligner(lam, path_model,supp_data=supp_data)
    }
    if(nrow(lam@peaks)==length(lam@clustering)){
      if(clustering){
        clustering <- finalClustering(lam,bw=lam@references@parameters$rt_dens,binSize=lam@references@parameters$dmz,bpp)
        lam <- clusterPeaktable(lam,clustering,bpp=bpp)
        lam@clustering <- clustering
        ###We just return the index
        return(lam)
      }else{
        lam@clustering <- 1:nrow(lam@peaks)
      }
      lam <- saveAligner(lam, path_model,supp_data=supp_data)
    }
    message("No more files to process.")
    return(lam)
  }

####We now perform the alignement by batch to speed the process

#' Align all the file form a directory to an exisiting LCMSAlignerModel eventually
#'
#' @param directory The directory to monitor
#' @param path_model The path to store the model
#' @param output The directory of output for for the ocnstructed block file
#' @param save_interval The number of aligned peak tbale between save
#' @param num_file If A new modle is created, how many files are used to find the references peaks
#' @param reset Shall the processing be restarted from the beginning eventually
#' @param span Smoothing parameters contorlling the loess
#' @param n_clusters Numbers of clusters evnetually tested.
#' @param ... Arguments which will be passed to LCMSAlignerModel
#'
#' @return The LCMSAlignerObjects
#' @export
#'
#' @examples
#' print("Examples to be put here")
LCMSAlignerModelFromDirectoryByBatch <-
  function(directory,
           path_model,
           output,
           save_interval = 5,
           size_batch = 50,
           num_file = 10,
           reset = FALSE,
           threshold = 5,
           span = 0.6,
           n_clusters = 15,
           supp_data=c("peakwidth","SN","right_on_left_assymetry"),
           ransac_l1=0.05,
           ransac_dist_threshold=0.1,
           maxAlign = 100,
           seed=512,
           max_cor = 0.1,
           clustering=TRUE,
           bpp=NULL,
           graphical=FALSE,
           ...) {
    message("Processing of directory ", directory, ".")

    if(length(save_interval)<size_batch) size_interval <- size_batch
    if(is.null(maxAlign)) maxAlign <- 100000
    if(is.null(bpp)) bpp <- bpparam()

    ###Two cases or the furnished vector is a set of chracter or the furnished data table is a set of character
    if(dir.exists(directory)){
      all_files <- list.files(directory, full.names = TRUE)
    }else{ ###It is a set fo files eventually.
      all_files <- directory
    }

    lam <- NULL
    ###If the file is null we intialize it.
    if (file.exists(path_model)) {
      lam <-
        LCMSAlignerModel(
          paths = path_model,
          save_interval = save_interval,
          output = output,
          n_clusters = n_clusters,
          supp_data=supp_data,
          graphical = graphical,
          bpp = bpp,
          ...
        )
      if (reset) {
        lam <- restart(lam,supp_data=supp_data)
      }
    } else{

      ###Change oa pradigm we tak the biggest ifle
      vsize <- file.info(all_files)$size
      ref_files <- all_files[order(vsize,decreasing=TRUE)[1:min(num_file, length(all_files))]]

      # ref_files <-
      #   sample(all_files, size = min(num_file, length(all_files)))


      lam <-
        LCMSAlignerModel(
          paths = ref_files,
          save_interval = save_interval,
          output = output,
          n_clusters = n_clusters,
          supp_data=supp_data,
          graphical = graphical,
          seed=seed,
          bpp = bpp,
          ...
        )
    }
    ###We check the file to be processed
    vpro <- isProcessed(lam, all_files)
    all_files <- all_files[!vpro]
    num_align <- 0
    while (length(all_files) != 0) {
      num_align <- num_align+1
      if(num_align>maxAlign) break
      all_files <- list.files(directory, full.names = TRUE)
      vpro <- isProcessed(lam, all_files)
      all_files <- all_files[!vpro]
      if (length(all_files) != 0) {
        message("Found ", length(all_files), " files to process.")
      }

      ###We select a batch of files to align.
      seq_batch <- seq(1,length(all_files),by = size_batch)
      nbatch <- length(seq_batch)
      if(seq_batch[length(seq_batch)] != length(all_files)){
        seq_batch <- c(seq_batch,length(all_files)+1)
      }

      for (nb in 1:nbatch) {
        ###We align each file individually
        lam <- alignPeaktableBatch(lam,
                              all_files[(seq_batch[nb]):(seq_batch[nb+1]-1)], path_aligner = path_model, span = span,
                              supp_data=supp_data,ransac_l1=ransac_l1,
                              ransac_dist_threshold = ransac_dist_threshold,
                              lim = max_cor, graphical = graphical, bpp = bpp)
      }
      ## we initilaize the clustering
      lam@clustering <- 1:nrow(lam@peaks)
    }
    ###We initialize the

    ##We do the ifnal clustering mode
    if(nrow(lam@data)!=0){
      lam <- saveAligner(lam, path_model,supp_data=supp_data)
    }
    if(nrow(lam@peaks)==length(lam@clustering)){
      if(clustering){
        clustering <- finalClustering(lam,bw=lam@references@parameters$rt_dens,
                                      binSize=lam@references@parameters$dmz,bpp)
        lam <- clusterPeaktable(lam,clustering,bpp=bpp)
        lam@clustering <- clustering
        ###We just return the index
        return(lam)
      }else{
        lam@clustering <- 1:nrow(lam@peaks)
      }
      lam <- saveAligner(lam, path_model,supp_data=supp_data)
    }
    message("No more files to process.")
    return(lam)
  }
