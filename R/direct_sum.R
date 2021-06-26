#' Direct Sum for Symmetric Matrices
#' 
#' Direct Sum for Symmetric Matrices
#' 
#' Computes the direct sum of the matrices in \code{x}
#' 
#' @param x list of matrices
#' 
#' @return Sparse Matrix
#' @import Matrix
#' 
#' @keywords internal
#' 
direct_sum <- function(x) {
  n <- length(x) 
  m <- sapply(x,nrow)
  m.cumulative <- apply(array(1:n),1,function(k){sum(m[1:k])})

  z <- expand.grid(col=1:m[1],row=1:m[1])
  out <- data.frame(row=z$row,col=z$col,value=as.vector(x[[1]]))
  if (n > 1) {
    for (i in 2:n) {
      z <- expand.grid(col=(1:m[i])+m.cumulative[i-1],row=(1:m[i])+m.cumulative[i-1])
      tmp <- data.frame(row=z$row,col=z$col,value=as.vector(x[[i]]))
      out <- rbind(out,tmp)
    }
  }
  out <- out[out$col >= out$row,]
  return(sparseMatrix(i=out[,1],j=out[,2],x=out[,3],symmetric=T))
}
