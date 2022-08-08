#' BLUP
#' 
#' BLUP
#' 
#' The argument \code{what} takes 5 possible values: "AV" (additive value), "BV" (breeding value), "GV" (genotypic value), "AM" (additive marker effect), and "DM" (dominance marker effect). "Values" refer to predictions for individuals, as opposed to markers. Predicted values include the average fixed effect of the environments, whereas predicted marker effects do not. Argument \code{index.weights} is a named vector (matching the names of the locations or traits), and the values are interpreted for standardized traits. 
#' 
#' @param data object of \code{\link{class_prep}} from \code{\link{blup_prep}}
#' @param geno object of \code{\link{class_geno}} from \code{\link{read_geno}}
#' @param what One of the following: AV, BV, GV, AM, DM. See Details.
#' @param index.weights named vector of index weights for the locations or traits
#' @param gwas.ncore Integer indicating number of cores to use for GWAS (default is 0 for no GWAS). 
#' 
#' @return Data frame of BLUPs
#' 
#' @import Matrix
#' @importFrom stats formula pnorm
#' @importFrom parallel makeCluster clusterExport parRapply stopCluster
#' @export

blup <- function(data,geno=NULL,what,index.weights=NULL,gwas.ncore=0L) {
  
  if (what %in% c("id","marker")) {
    stop("The options for 'what' have changed. See Details.")
  }
  what <- toupper(what)
  stopifnot(what %in% c("AV","BV","GV","AM","DM"))
  
  if (nrow(data@add) > 0 & is.null(geno))
    stop("Missing geno argument")
  
  if (nrow(data@add)==0 & (!is.null(geno) | what %in% c("AM","DM","AV","BV")))
    stop("Marker data was not used in blup_prep")
  
  if (what=="DM" & (!(class(geno)=="class_genoD") | nrow(data@dom)==0))
    stop("Dominance model is required in geno and data")
  
  gamma <- 0
  if (!is.null(geno)) {
    if (what=="GV")  
      gamma <- 1
    if (what=="BV" & nrow(data@dom) > 0)
      gamma <- (geno@ploidy/2 - 1)/(geno@ploidy - 1)
    
    if (nrow(data@dom) > 0) {
      index.scale <- sqrt(diag(as.matrix(data@add) + gamma^2*diag(as.matrix(data@dom))))
    } else {
      index.scale <- sqrt(diag(as.matrix(data@add) + gamma^2*diag(as.matrix(data@g.iid))))
    }
  } else {
    index.scale <- sqrt(diag(as.matrix(data@g.iid)))
  }
  
  n.id <- length(data@id)
  nlt <- length(index.scale) #number of loc/traits
  
  if (nlt > 1) {
    if (is.null(index.weights)) {
      index.coeff <- 1/index.scale
    } else {
      ix <- match(names(index.scale),names(index.weights))
      if (any(is.na(ix))) {stop("Check names of loc or traits in index.weights")}
      index.weights <- index.weights[ix]
      index.coeff <- index.weights/index.scale
    }
    index.coeff <- index.coeff/sqrt(sum(index.coeff^2))  
    fix.value <- sum(index.coeff*data@avg.env[names(index.coeff)])
  } else {
    index.coeff <- 1
    fix.value <- data@avg.env
  }
  
  n.mark <- length(data@fixed.marker)/nlt
  if (n.mark > 0) {
    if (nlt > 1) {
      fixed.markers <- unique(sapply(strsplit(names(data@fixed.marker),split=":",fixed=T),"[[",1))
      beta <- matrix(data@fixed.marker,ncol=1)
      tmp <- kronecker(geno@coeff[,fixed.markers,drop=FALSE],matrix(index.coeff,nrow=1)) %*% beta
      fix.value <- fix.value + as.numeric(tmp)
      
    } else {
      fixed.markers <- names(data@fixed.marker)
      beta <- matrix(data@fixed.marker,ncol=1)
      fix.value <- fix.value + as.numeric(geno@coeff[,fixed.markers,drop=FALSE] %*% beta)
    } 
  }
  
  if (length(data@heterosis) > 0 & what %in% c("BV","GV")) {
    fix.value <- fix.value + gamma * geno@Fg[data@id] * sum(index.coeff*data@heterosis)
  }
  
  if (what %in% c("AV","BV","GV")) {
    if (gwas.ncore > 0)
      stop("GWAS option requires either AM or DM for 'what' ")
    
    M <- kronecker(Diagonal(n=n.id),Matrix(index.coeff,nrow=1))
    if (!is.null(geno)) {
      M <- cbind(M,gamma*M)
    }
    out <- data.frame(id=data@id,value=as.numeric(M%*%Matrix(data@random,ncol=1)) + fix.value)
    numer <- diag(M%*% tcrossprod(data@var.uhat,M))
    denom <- diag(M%*% tcrossprod(data@var.u,M))
    out$r2 <- numer/denom
    attr(out,"what") <- what
    rownames(out) <- NULL
    return(out)
  }
  
  #marker effects
  n.random <- length(data@random)
  if (!is.null(geno)) 
    n.random <- n.random/2
  
  if (what=="AM") {
    G <- kron(geno@eigen.G,1)
    M <- kronecker(crossprod(geno@coeff,crossprod(G$inv)),matrix(index.coeff,nrow=1))/geno@scale
    effect <- as.numeric(M %*% matrix(data@random[1:n.random],ncol=1))
    V <- data@var.uhat[1:n.random,1:n.random] #used for GWAS
  } else {
    D <- kron(geno@eigen.D,1)
    M <- kronecker(crossprod(geno@coeff.D,crossprod(D$inv)),matrix(index.coeff,nrow=1))/geno@scale.D
    effect <- as.numeric(M %*% matrix(data@random[n.random + 1:n.random],ncol=1))
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
