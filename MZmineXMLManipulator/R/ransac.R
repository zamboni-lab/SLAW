
fitLoess <- function(x,y,span=0.4){
  rloess <- loess(y~x,span = span)
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
  mean(calculateDistancePred(y,predy))**2+alpha*mean(abs(predy*vpen))
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
  legend("topleft",legend=c("Model seeds","Inliers","Miss detected"),
          col=c("orange","blue","black"),pch=1)
  return(model)
}


ransacLoess <- function(tdata,max_iter=100,fitting_point=floor(nrow(tdata)*0.1),
                        dist_threshold=5,min_inliers=50,graphical=FALSE,span=0.3,L1=0.05){
  inliers <-  NULL
  nits <-  0
  best_model <-  NULL
  best_error <- 100000
  ndata <- nrow(tdata)
  allinliers <- NULL
  binliers <- NULL
  bpoints <- NULL
  sel_points <- NULL

  seqx <- seq(min(tdata$x),max(tdata$x),length=15)

  while(nits < max_iter){
    ###We sample a ramdon numbr of data points
    sel_points <- sample.int(ndata,fitting_point)

    ###We fit the model
    current_model <- fitLoess(tdata$x[sel_points],tdata$y[sel_points],span=span)
    my <- max(tdata$y[sel_points])
    my <- mean(abs(tdata$y[sel_points]))*4
    ##At each point we check if there is no abbherrant point
    vpred <- tryCatch(predict(current_model,seqx),error=function(e) return(NA))
    if(length(vpred)==1&&is.na(vpred)) next
    vpred <- vpred[!is.na(vpred)]
    if(any(vpred>my)){
      next
    }

    ###We calculate the distance for all the points
    current_distance <- tryCatch(calculateDistance(current_model,tdata$x,tdata$y),error=function(e){return(NA)})
    if(length(current_distance)==1&&is.na(current_distance)) next

    ###we calculate the best residuals distance
    inliers <- which(current_distance<dist_threshold)

    ###
    inliers <- setdiff(inliers,sel_points)

    if(length(inliers)>min_inliers){
      ##We first hceck the max amplitude of th emodel on 10 points


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
      }
    }
    nits <- nits+1
  }
  if(graphical){
    best_model <- tryCatch(plotModel(best_model,tdata,binliers,bpoints),error=function(e){return(NA)})
    #tryCatch(plotModel(best_model,tdata,binliers,bpoints),error=function(e){return(NA)})
  }
  return(best_model)
}

####
