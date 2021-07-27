#' Create Manhattan plot
#' 
#' Create Manhattan plot
#' 
#' Assumes position in bp
#' 
#' @param data data frame with columns for marker, chrom, position, and gwas.score
#' @param chrom optional, to plot only one chromosome
#' @param thresh optional, to include horizontal line at discovery threshold
#' @param rotate TRUE/FALSE whether to rotate x-axis labels to be perpendicular 
#' 
#' @return ggplot2 object
#' 
#' @export
#' @import ggplot2
#' @importFrom rlang .data

manhattan_plot <- function(data,chrom=NULL,thresh=NULL,rotate=FALSE) {
		
  stopifnot(c("chrom","position","gwas.score") %in% colnames(data))
  if (is.null(chrom)) {
    x <- get_x(data[,c("chrom","position")])
    ix <- 1:nrow(data)
  } else {
    stopifnot(chrom %in% unique(data$chrom))
    ix <- which(as.character(data$chrom)==chrom)
    x <- data$position[ix]/1e6
  }
  plot.data <- data.frame(x=x,y=data$gwas.score[ix],
                          color=factor(ifelse(as.integer(factor(data$chrom[ix]))%%2==1,1,0)))
  if (rotate) {
    angle <- 90
  } else {
    angle <- 0
  }
  p <- ggplot(data=plot.data,
              aes(x=.data$x,y=.data$y,colour=.data$color)) +
    ylab(expression(paste("-log"[10],"(p)"))) + guides(colour="none") + 
    theme_bw() + theme(text = element_text(size=15),axis.text.x=element_text(angle=angle,vjust=0.5)) + geom_point() 
  
  if (is.null(chrom)) {
    allchr <- unique(data$chrom)
    breaks <- (tapply(x,data$chrom,max) + tapply(x,data$chrom,min))/2
    p <- p + scale_x_continuous(name="Chromosome",breaks=breaks,labels=allchr) +
     scale_colour_manual(values=c("#21908c","#440154"))
  } else {
    p <- p + scale_x_continuous(name="Position (Mb)") + scale_colour_manual(values="#440154")
  }
    
  if (!is.null(thresh)) {
    p <- p + geom_hline(yintercept = thresh,linetype=2,colour="grey50")
  }
	return(p)	
}
