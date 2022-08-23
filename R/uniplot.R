#' Uniplot for multi-location models
#'
#' Displays scaled loadings of the FA2 model
#' 
#' The squared radius for each point is the proportion of genetic variance explained by the latent factors. For points on the unit circle, the cosine of the subtended angle equals the correlation.
#'
#' @param loadings scaled factor loadings, from \code{\link{Stage2}}. 
#' @param nudge distance to nudge labels
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @return ggplot2 object
#' @export

uniplot <- function(loadings,nudge=0.1) {
  x <- seq(-1,1,by=0.01)
  y1 <- sqrt(1-x^2)
  y2 <- -y1
  plot.data <- data.frame(name=rownames(loadings),x=loadings[,1],y=loadings[,2])
  p <- ggplot(data=data.frame(x,y1,y2)) + geom_line(mapping=aes(x=.data$x,y=.data$y1),colour="red") + geom_line(mapping=aes(x=.data$x,y=.data$y2),colour="red") + theme_bw() + coord_fixed(ratio=1) +
    geom_point(mapping=aes(x=.data$x,y=.data$y),data=plot.data) + xlab("") + ylab("") + 
    theme(axis.text = element_blank(),axis.ticks = element_blank()) 
  p + geom_text_repel(aes(x=.data$x,y=.data$y,label=.data$name),force=2,
                      nudge_x=nudge,plot.data,segment.colour="blue")
}