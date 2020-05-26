
fitLoess <- function(x,y,weights=NULL,span=0.4){
  if(is.null(weights)) weights <- rep(1,length(y))
  rloess <- suppressWarnings(loess(y~x,weights = weights,span = span))
  ###we return the model
  return(rloess)
}


calculateDistance <- function(model,x,y){
  abs(y-predict(model,newdata=x))
}

calculateDistancePred <- function(y,predy){
  abs(y-predy)
}

scoreModel <- function(model,x,y,alpha=0.05){
  ###We add a specfic term to reduce th eimpact of early deviation.
  predy <- predict(model,newdata=x)
  vpen <- (1-0.8*x/max(x))
  mean(calculateDistancePred(y,predy)**2)+alpha*mean(abs(predy*vpen))
}


plotModel <-  function(model,tdata,inliers,inpoints){
  col <- rep("black",nrow(tdata))
  col[inliers] <- "blue"
  col[inpoints] <-  "orange"
  ssx <- seq(min(tdata$x),max(tdata$x),length=300)
  allp <- predict(model,ssx)
  ssx <- ssx[!is.na(allp)]
  allp <- allp[!is.na(allp)]

  ylim <- max(c(allp,tdata$y))
  ylim <- c(-ylim,ylim)
  plot(tdata$x,tdata$y,col=col,type="p",xlab="RT",ylab=c("RT deviation"),ylim=ylim)
  # browser()
  ##We add the model

  lines(ssx,allp,col="red")
  # legend("topleft",legend=c("Model seeds","Inliers","Miss detected"),
  #         col=c("orange","blue","black"),pch=1)
  return(model)
}



loRansacLoss <- function(tdata,max_iter=100,fitting_point=floor(nrow(tdata)*0.1),
                         dist_threshold=5,min_inliers=50,graphical=FALSE,span=0.3,L1=0.05,
                         weight=NULL,max_iter_sample=3,max_optim=5,sub_optim=5){
  
  if(is.null(weight)) weight <- rep(1,nrow(tdata))
  inliers <-  NULL
  nits <-  0
  best_model <-  NULL
  # best_error <- 100000
  ndata <- nrow(tdata)
  allinliers <- NULL
  binliers <- NULL
  bpoints <- NULL
  sel_points <- NULL
  iter_sample <- 0
  
  seqx <- seq(min(tdata$x),max(tdata$x),length=15)
  pleft <- c(seqx[1],median(tdata$y[tdata$x<=seqx[1]]))
  pright <- c(seqx[15],median(tdata$y[tdata$x>=seqx[14]]))
  
  valid_model <- 0
  inlierscount <- rep(1,nrow(tdata))
  # browser()
  while(iter_sample<max_iter_sample){
    probas <- weight/sum(weight)
    mprob <- min(probas)
    while(nits < max_iter){
      ###We sample a ramdon numbr of data points
      sel_points <- sample.int(ndata,fitting_point,prob = probas)
      nits <- nits+1
      ###We fit the model
      current_model <- fitLoess(tdata$x[sel_points],tdata$y[sel_points],weights=probas[sel_points],span=span)
      my <- max(tdata$y[sel_points])
      my <- mean(abs(tdata$y[sel_points]))*4
      ##At each point we check if there is no abbherrant point
      vpred <- tryCatch(predict(current_model,seqx),error=function(e) return(NA))
      if(length(vpred)==1&&is.na(vpred)) next
      vpred <- vpred[!is.na(vpred)]
      if(any(vpred>my)) next
      
      ###We calculate the distance for all the points
      current_distance <- tryCatch(calculateDistance(current_model,tdata$x,tdata$y),error=function(e){return(NA)})
      if(length(current_distance)==1&&is.na(current_distance)) next
      
      ###we calculate the best residuals distance
      inliers <- which(current_distance<dist_threshold)
      o_inliers <- inliers
      inliers <- setdiff(inliers,sel_points)
      
      if(length(inliers)>min_inliers){
        
        ##We perform inner optimization of the parameters
        num_optim <- 0
        best_inliers_optim <- inliers
        best_model_optim <- current_model
        best_point_optim <- sel_points
        o_inliers_optim <- o_inliers
        while(num_optim<max_optim){
          num_optim <- num_optim+1
          fitting_optim <- max(min_inliers,length(inliers)/2)
          sel_points_optim <- sample(o_inliers,size = fitting_optim,prob = probas[o_inliers])
          ###We add the minimum points

          ###We add the rest of the data.
          ###We fit the model
          optim_model <- fitLoess(c(tdata$x[sel_points_optim],pleft[1],pright[1]),c(tdata$y[sel_points_optim],pleft[2],pright[2]),
                                  weights=c(mprob,probas[sel_points_optim],mprob),span=span)
          my <- max(tdata$y[sel_points_optim])
          my <- mean(abs(tdata$y[sel_points_optim]))*4
          ##At each point we check if there is no abbherrant point
          vpred <- tryCatch(predict(optim_model,seqx),error=function(e) return(NA))
          if(length(vpred)==1&&is.na(vpred)) next
          vpred <- vpred[!is.na(vpred)]
          if(any(vpred>my)) next
          
          ###We calculate the distance for all the points
          current_distance_optim <- tryCatch(calculateDistance(optim_model,tdata$x,tdata$y),error=function(e){return(NA)})
          if(length(current_distance_optim)==1&&is.na(current_distance_optim)) next
          
          ###we calculate the optimized inliers count
          inliers_optim <- which(current_distance_optim<dist_threshold)
          o_inliers_optim <- inliers_optim
          inliers_optim <- setdiff(inliers_optim,sel_points_optim)
          if(length(inliers_optim)>length(best_inliers_optim)){
            # message("OPTIM")
            best_inliers_optim <- inliers_optim
            best_model_optim <- optim_model
            best_point_optim <- sel_points_optim
          }
        }
        inlierscount[best_inliers_optim] <- inlierscount[best_inliers_optim]+1
        valid_model <- valid_model+1
        best_model <- best_model_optim
        binliers <- o_inliers_optim
        bpoints <- best_point_optim
        if(graphical){
          plotModel(best_model,tdata,binliers,bpoints)
        }
      }
    }
    weight <- inlierscount/sum(inlierscount)
    iter_sample <- iter_sample+1
  }
  if(graphical){
    best_model <- tryCatch(plotModel(best_model,tdata,binliers,bpoints),error=function(e){return(NA)})
    #tryCatch(plotModel(best_model,tdata,binliers,bpoints),error=function(e){return(NA)})
  }
  # message("valid_model:",valid_model)
  return(best_model)
}


ransacLoess <- function(tdata,max_iter=100,fitting_point=floor(nrow(tdata)*0.1),
                        dist_threshold=5,min_inliers=50,graphical=FALSE,span=0.3,L1=0.05,weight=NULL,max_iter_sample=3){
  
  if(is.null(weight)) weight <- rep(1,nrow(tdata))
  inliers <-  NULL
  nits <-  0
  best_model <-  NULL
  best_error <- 100000
  ndata <- nrow(tdata)
  allinliers <- NULL
  binliers <- NULL
  bpoints <- NULL
  sel_points <- NULL
  iter_sample <- 0

  seqx <- seq(min(tdata$x),max(tdata$x),length=15)
  valid_model <- 0
  inlierscount <- rep(1,nrow(tdata))
  # browser()
  while(iter_sample<max_iter_sample){
    probas <- weight/sum(weight)
    # boxplot(weight)
    while(nits < max_iter){
      ###We sample a ramdon numbr of data points
      sel_points <- sample.int(ndata,fitting_point,prob = probas)
      nits <- nits+1
      ###We fit the model
      current_model <- fitLoess(tdata$x[sel_points],tdata$y[sel_points],weights=probas[sel_points],span=span)
      my <- max(tdata$y[sel_points])
      my <- mean(abs(tdata$y[sel_points]))*4
      ##At each point we check if there is no abbherrant point
      vpred <- tryCatch(predict(current_model,seqx),error=function(e) return(NA))
      if(length(vpred)==1&&is.na(vpred)) next
      vpred <- vpred[!is.na(vpred)]
      if(any(vpred>my)) next
  
      ###We calculate the distance for all the points
      current_distance <- tryCatch(calculateDistance(current_model,tdata$x,tdata$y),error=function(e){return(NA)})
      if(length(current_distance)==1&&is.na(current_distance)) next
  
      ###we calculate the best residuals distance
      inliers <- which(current_distance<dist_threshold)
      o_inliers <- inliers
      inliers <- setdiff(inliers,sel_points)
  
      if(length(inliers)>min_inliers){
        inlierscount[o_inliers] <- inlierscount[o_inliers]+1
        valid_model <- valid_model+1
  
        ###we check if the best error is smaller
        allinliers <- c(inliers,sel_points)
        current_error <- tryCatch(scoreModel(current_model,tdata$x[allinliers],tdata$y[allinliers],alpha=L1),error=function(e){
          return(NA)
        })
        if(length(current_error)==1&&is.na(current_error)) next
        if(current_error<best_error){
          best_error <- current_error
          best_model <- current_model
          binliers <- allinliers
          bpoints <- sel_points
          if(graphical){
          plotModel(best_model,tdata,binliers,bpoints)
          }
        }
      }
    }
    weight <- inlierscount/sum(inlierscount)
    iter_sample <- iter_sample+1
  }
  if(graphical){
    best_model <- tryCatch(plotModel(best_model,tdata,binliers,bpoints),error=function(e){return(NA)})
    #tryCatch(plotModel(best_model,tdata,binliers,bpoints),error=function(e){return(NA)})
  }
  message("valid_model:",valid_model)
  return(best_model)
}

####
