#' Remove spatial trend 
#' 
#' Removes spatial trend to prepare for multi-trait Stage 1
#' 
#' SpATS used to remove 2D spatial trend for each field experiment, labeled with column 'expt' in the input file. Genotype labels are in column 'id' and modeled as i.i.d random effect. Argument \code{effects} is used to specify other i.i.d. effects besides genotype and has three columns: name, fixed, factor. The "name" column is a string that must match a column in the input file. The fixed column is a logical variable to indicate whether the effect is fixed (TRUE) or random (FALSE). The factor column is a logical variable to indicate whether the effect is a factor (TRUE) or numeric (FALSE). 
#' 
#' Argument \code{traits} is a character vector of trait names. Single-trait analyses are performed for each trait, and the results are combined in the output.
#' 
#' @param filename Name of CSV input file
#' @param traits trait names 
#' @param spline vector of variable names for 2D spline with SpATS
#' @param effects data frame specifying other effects in the model (see Details)
#' 
#' @return Data frame with adjusted phenotypes
#' 
#' @importFrom utils read.csv
#' @import SpATS
#' @export

remove_spatialtrend <- function(filename, traits, spline, effects=NULL) {
  
  data <- read.csv(file=filename,check.names=F)
  stopifnot(traits %in% colnames(data))
  stopifnot(c("id","expt") %in% colnames(data))
  n.trait <- length(traits)
  
  if (!is.null(effects)) {
    stopifnot(effects$names %in% colnames(data))
  }
  
  data$expt <- as.character(data$expt)
  data$id <- as.character(data$id)
  expts <- unique(data$expt)
  n.expt <- length(expts)

  model <- "SpATS(data=data1,response='y',genotype='id',fixed=~FIX,random=~RANDOM,spatial=~SAP(spline.x,spline.y),genotype.as.random=TRUE,control=list(monitoring=0))"
  
  model <- sub("spline.x",spline[1],model,fixed=T)
  model <- sub("spline.y",spline[2],model,fixed=T)
  
  if (!is.null(effects)) {
    k <- which(effects$fixed)
    if (length(k) > 0) {
      FIX <- paste(effects$name[k],collapse="+")
      model <- sub("FIX",FIX,model,fixed=T)
    } else {
      model <- sub("fixed=~FIX,","",model,fixed=T)
    }
    k <- which(effects$random)
    if (length(k) > 0) {
      RANDOM <- paste(effects$name[k],collapse="+")
      model <- sub("RANDOM",RANDOM,model,fixed=T)
    } else {
      model <- sub("random=~RANDOM,","",model,fixed=T)
    }
  } else {
    model <- sub("fixed=~FIX,","",model,fixed=T)
    model <- sub("random=~RANDOM,","",model,fixed=T)
  }
  
  adjusted <- NULL
  for (i in 1:n.expt) {
    ix <- which(data$expt==expts[i])
    data1 <- data[ix,]
    nx <- length(unique(data1[,spline[1]]))
    ny <- length(unique(data1[,spline[2]]))
    
    if (!is.null(effects)) {
      k <- which(effects$factor)
      if (length(k)>0) {
        for (q in 1:k)
          eval(parse(text="data1[,effects$name[q]] <- factor(as.character(data1[,effects$name[q]]))"))
      }
    }
    
    data1$id <- factor(data1$id)
    for (j in 1:n.trait) {
      data1$y <- data1[,traits[j]]
      ans <- eval(parse(text=model))
      st <- obtain.spatialtrend(ans,grid=c(nx,ny))
      st1 <- data.frame(expand.grid(x=st$col.p, y=st$row.p),
                        z=as.numeric(t(st$fit)))
      colnames(st1) <- c(spline,"spatialtrend")
      data2 <- merge(data1,st1,by=spline)
      data2[,traits[j]] <- data2$y - data2$spatialtrend
      if (j==1) {
        out <- data2[,-match(c("y","spatialtrend"),colnames(data2))]
      } else {
        out <- merge(out[,-match(traits[j],colnames(out))],
                     data2[,-match(c("y","spatialtrend",traits[-j]),colnames(data2))])
      }
    }
    adjusted <- rbind(adjusted,out)
  }
  return(adjusted)
}

