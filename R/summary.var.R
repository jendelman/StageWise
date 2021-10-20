#' Summarize variances and correlations
#'
#' Summarize variances and correlations
#' 
#' When reporting the partitioning of variance, the variance component for the additive effect is multiplied by the mean diagonal of the G matrix. When a G matrix has been included, the correlation matrix shown is for the additive value (plus any markers).
#' 
#' @param object object of \code{\link{class_var}}
#'
#' @return For univariate analysis, a matrix with the proportion of variance. For multi-location or multi-trait analysis, a list containing the proportion of variance and the correlation matrix.
#'
#' @include class_var.R
#' @name summary
#' @exportMethod summary
NULL

setGeneric("summary")
setMethod("summary",c(object="class_var"),
          definition=function(object){
            cor.mat <- NULL
            if (ncol(object@resid)>1) {
              traits <- rownames(object@resid)
              #multiple traits
              if (is.na(object@meanG)) {
                cor.mat <- round(cov_to_cor(object@g.resid),2)
                vtab <- data.frame(matrix(c(diag(object@g.resid),diag(object@resid)),nrow=2,byrow = T))
                colnames(vtab) <- traits
                if (any(!is.na(object@meanOmega))) {
                  vtab <- rbind(vtab,object@meanOmega)
                  rownames(vtab) <- c("genotype","g x env","residual")
                } else {
                  rownames(vtab) <- c("genotype","residual")
                }
              } else {
                add.cov <- object@add
                if (nrow(object@fixed.marker.cov)>0)
                  add.cov <- add.cov + object@fixed.marker.cov
                cor.mat <- round(cov_to_cor(add.cov),2)
                vtab <- data.frame(matrix(c(diag(object@add)*object@meanG,diag(object@g.resid),diag(object@resid)),nrow=3,byrow = T))
                colnames(vtab) <- traits
                if (any(!is.na(object@meanOmega))) {
                  vtab <- rbind(vtab,object@meanOmega)
                  rownames(vtab) <- c("additive","g.resid","g x env","residual")
                } else {
                  rownames(vtab) <- c("additive","g.resid","residual")
                }
              }
              
              if (nrow(object@fixed.marker.var)>0) {
                  vtab <- rbind(object@fixed.marker.var,vtab)
              }
              vtab <- as.matrix(vtab)
              pvar <- vtab %*% diag(1/apply(vtab,2,sum))
              colnames(pvar) <- traits
              return(list(round(pvar,3),cor.mat))
            }
            
            #single trait
            if (!is.na(object@meanG)) {
              if (nrow(object@add) > 1) {
                #multiple locations
                add.cov <- object@add
                if (nrow(object@fixed.marker.cov)>0)
                  add.cov <- add.cov + object@fixed.marker.cov
                cor.mat <- round(cov_to_cor(add.cov),2)
                tmp <- partition(object@add)
                Va <- tmp[1]*object@meanG
                VaL <- tmp[2]*object@meanG
                
                Vr <- mean(object@g.resid[,1])
                V1 <- matrix(c(Va,VaL,Vr),ncol=1)
                rownames(V1) <- c("additive","add x loc","g.resid")
              } else {
                #one location
                Va <- as.numeric(object@add)*object@meanG
                Vr <- as.numeric(object@g.resid)
                
                V1 <- matrix(c(Va,Vr),ncol=1)
                rownames(V1) <- c("additive","g.resid")
              }
            } else {
              if (nrow(object@g.resid) > 1) {
                cor.mat <- round(cov_to_cor(object@g.resid),2)
                Vg <- mean(object@g.resid[upper.tri(object@g.resid,diag=F)])
                VgL <- mean(diag(object@g.resid)) - Vg
                
                V1 <- matrix(c(Vg,VgL),ncol=1)
                rownames(V1) <- c("genotype","g x loc")
              } else {
                Vg <- as.numeric(object@g.resid)
                V1 <- matrix(c(Vg),ncol=1)
                rownames(V1) <- c("genotype")
              }
            }
            
            if (!is.na(object@meanOmega)) {
              V2 <- matrix(c(mean(object@resid),object@meanOmega),ncol=1)
              rownames(V2) <- c("g x env","residual")
            } else {
              V2 <- matrix(mean(object@resid),ncol=1)
              rownames(V2) <- "residual"
            }
            vtab <- rbind(V1,V2)
            
            if (nrow(object@fixed.marker.var)>0) {
              if (ncol(object@fixed.marker.var)>1) {
                markers <- rownames(object@fixed.marker.var)
                x <- matrix(t(object@fixed.marker.var),ncol=1)
                tmp <- expand.grid(factor(markers,levels = markers,ordered=T),c(""," x loc"))
                rownames(x) <- apply(tmp[order(tmp$Var1),],1,paste,collapse="")
              } else {
                x <- object@fixed.marker.var
              }
              vtab <- rbind(x,vtab)
            }
            
            pvar <- round(vtab / sum(vtab),3)
            colnames(pvar) <- "Var. Proportion"
            if (!is.null(cor.mat)) {
              return(list(pvar,as.matrix(cor.mat)))
            } else {
              return(pvar)
            }
          })