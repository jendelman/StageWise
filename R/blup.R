#' BLUP
#' 
#' BLUP
#' 
#' The argument \code{what} takes 5 possible values: "AV" (additive value), "BV" (breeding value), "GV" (genotypic value), "AM" (additive marker effect), and "DM" (dominance marker effect). "Values" refer to predictions for individuals, as opposed to markers. Predicted values include the average fixed effect of the environments, whereas predicted marker effects do not. Argument \code{index.coeff} is a named vector (matching the names of the locations or traits), and the values are interpreted for standardized traits. 
#' 
#' When multiple objects of \code{\link{class_prep}} are used for \code{data}, they must be based on the same marker data and genetic model. Also, reliabilities are not computed.
#' 
#' @param data one object, or list of objects, of \code{\link{class_prep}} from \code{\link{blup_prep}}
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param what One of the following: AV, BV, GV, AM, DM. See Details.
#' @param index.coeff named vector of index coefficients for the locations or traits
#' @param gwas.ncore Integer indicating number of cores to use for GWAS (default is 0 for no GWAS). 
#' 
#' @return Data frame of BLUPs
#' 
#' @import Matrix
#' @importFrom stats formula pnorm
#' @importFrom parallel makeCluster clusterExport parRapply stopCluster
#' @export

blup <- function(data, geno=NULL, what, index.coeff=NULL, gwas.ncore=0L) {
  
  if (what %in% c("id","marker")) {
    stop("The options for 'what' have changed. See Details.")
  }
  what <- toupper(what)
  stopifnot(substr(what,1,2) %in% c("AV","BV","GV","AM","DM"))
  
  if (is(data,"list")) {
    stopifnot(sapply(data,class)=="class_prep")
    multi.prep <- TRUE
  } else {
    stopifnot(is(data,"class_prep"))
    multi.prep <- FALSE
    data2 <- data
    data <- list(data)
  }
  
  if (data[[1]]@model > 0L & is.null(geno))
    stop("Missing geno argument")
  
  if (data[[1]]@model==0L & (!is.null(geno) | what %in% c("AM","DM","AV","BV")))
    stop("Marker data was not used in blup_prep")
  
  if (substr(what,1,2)=="DM" & (!(is(geno,"class_genoD")) | data[[1]]@model < 3L))
    stop("Dominance model is required in geno and data")
  
  n.id <- length(data[[1]]@id)
  gamma <- 0
  if (!is.null(geno)) {
    if (substr(what,1,2)=="GV")  
      gamma <- 1
    if (substr(what,1,2)=="BV" & data[[1]]@model==3L)
      gamma <- (geno@ploidy/2 - 1)/(geno@ploidy - 1)
    
    if (length(gamma) > 1) {
      tmp <- tapply(gamma,attr(geno@ploidy,"pop"),mean)
      gamma <- sum(tmp[names(index.coeff)]*index.coeff)/sum(index.coeff)
    }
  }
  
  if (data[[1]]@model < 2L) {
    index.scale <- sqrt(unlist(sapply(data,function(z){diag(as.matrix(z@geno1.var))})))
  } else {
    index.scale <- sqrt(unlist(sapply(data,function(z){diag(as.matrix(z@geno1.var))})) + 
                          gamma^2*unlist(sapply(data,function(z){diag(as.matrix(z@geno2.var))})))
  } 
  
  if (!multi.prep) {
    data <- data2
    names(index.scale) <- rownames(data@geno1.var)
    rm("data2")
  } else {
    nt <- sapply(data,function(z){nrow(z@geno1.var)})
    ix <- which(nt==1)
    if (length(ix) > 0) {
      for (k in ix)
        rownames(data[[k]]@geno1.var) <- names(data)[k]
    }
    names(index.scale) <- unlist(sapply(data,function(z){rownames(z@geno1.var)}))
  }
  
  
  nlt <- length(index.scale) #number of loc/traits
  
  if (nlt > 1) {
    traits <- names(index.scale)
    ix <- match(traits,names(index.coeff))
    if (any(is.na(ix))) {stop("Check names of loc/traits in index.coeff")}
    index.coeff <- index.coeff[ix]
    
    if (nchar(what)==2) {
      index.coeff <- index.coeff/as.numeric(index.scale)
      index.coeff <- index.coeff/sqrt(sum(index.coeff^2))
    } else {
      #do not rescale
      what <- substr(what,1,2)
    }
    names(index.coeff) <- traits
  } else {
    if (nchar(what)==2) {
      index.coeff <- 1
    } else {
      what <- substr(what,1,2)
    }
  }
  
  if (multi.prep) {
    for (i in 1:length(data)) {
      ix <- match(rownames(data[[i]]@geno1.var),names(index.coeff))
      tmp <- blup(data[[i]],geno=geno,what=paste0(what,"*"),index.coeff=index.coeff[ix])
      if (i==1) {
        if (what %in% c("AM","DM")) {
          out <- tmp$effect
        } else {
          out <- tmp[,1:2]
        }
      } else {
        if (what %in% c("AM","DM")) {
          out <- out + tmp$effect
        } else {
          out <- merge(out,tmp[,1:2],by="id")
        }
      }
    }
    if (what %in% c("AM","DM")) {
      tmp$effect <- out
      return(tmp)
    } else {
      return(data.frame(id=out$id, value=apply(out[,-1],1,sum)))
    }
  }
  
  if (nlt > 1) {
    fix.value <- sum(index.coeff*data@avg.env[names(index.coeff)])
  } else {
    fix.value <- data@avg.env*index.coeff
  }
  
  n.mark <- length(data@fixed.marker)/nlt
  if (n.mark > 0) {
    if (nlt > 1) {
      fixed.markers <- unique(sapply(strsplit(names(data@fixed.marker),split=":",fixed=T),"[[",1))
    } else {
      fixed.markers <- names(data@fixed.marker)
    }
    beta <- matrix(data@fixed.marker,ncol=1)
    tmp <- kronecker(geno@coeff[,fixed.markers,drop=FALSE],matrix(index.coeff,nrow=1)) %*% beta
    fix.value <- fix.value + as.numeric(tmp)
  }
  
  if (length(data@heterosis) > 0 & what %in% c("BV","GV")) {
    fix.value <- fix.value + gamma * (-geno@Fg[data@id]) * sum(index.coeff*data@heterosis)
  }
  
  if (what %in% c("AV","BV","GV")) {
    if (gwas.ncore > 0)
      stop("GWAS option requires either AM or DM for 'what' ")
    
    M <- kronecker(Diagonal(n=n.id),Matrix(index.coeff,nrow=1))
    if (data@model > 1L) 
      M <- cbind(M,gamma*M)
    
    out <- data.frame(id=data@id,value=as.numeric(M%*%Matrix(data@random,ncol=1)) + fix.value)
    numer <- diag(M%*% tcrossprod(data@var.uhat,M))
    denom <- diag(M%*% tcrossprod(data@var.u,M))
    out$r2 <- numer/denom
    attr(out,"what") <- what
    rownames(out) <- NULL
    return(out)
  }
  
  #marker effects
  if (data@model > 1L) {
    n.random <- length(data@random)/2
  } else {
    n.random <- length(data@random)
  }
  
  if (what=="AM") {
    G <- kron(geno@eigen.G,1)
    M <- kronecker(crossprod((geno@coeff/geno@scale),crossprod(G$inv)),matrix(index.coeff,nrow=1))
    effect <- as.numeric(M %*% matrix(data@random[1:n.random],ncol=1))
    V <- data@var.uhat[1:n.random,1:n.random] #used for GWAS
  } else {
    D <- kron(geno@eigen.D,1)
    M <- kronecker(crossprod((geno@coeff.D/geno@scale.D),crossprod(D$inv)),matrix(index.coeff,nrow=1))
    dhat <- -kronecker(matrix(geno@Fg[data@id],ncol=1),matrix(data@heterosis,ncol=1)) + 
      matrix(data@random[n.random + 1:n.random],ncol=1)
    effect <- as.numeric(M %*% dhat)
    V <- data@var.uhat[n.random+1:n.random,n.random+1:n.random] #used for GWAS
  }
  
  if (nrow(geno@map)==0) {
      out <- data.frame(marker=colnames(geno@coeff), effect=effect)
  } else {
      out <- data.frame(geno@map, effect=effect)
  }
  
  f.se <- function(x,V) {sqrt(as.numeric(crossprod(x,V%*%x)))}
  if (gwas.ncore > 0) {

    if (gwas.ncore == 1) {
      se <- apply(M,1,f.se,V=V)
      std.effect <- out$effect/se
    } else {
      cl <- makeCluster(gwas.ncore)
      clusterExport(cl=cl,varlist=NULL)
      se <- parRapply(cl=cl,x=M,f.se,V=V)
      stopCluster(cl)
      std.effect <- out$effect/sapply(se,as.numeric)
    } 
    out$score <- -log10(pnorm(q=abs(std.effect),lower.tail=FALSE)*2)
  }
  
  if (n.mark > 0) {
    k <- match(fixed.markers,out$marker)
    out$effect[k] <- out$effect[k] + as.numeric(beta)
  }
  attr(out,"what") <- what
  rownames(out) <- NULL
  return(out)
}
