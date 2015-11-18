#########
# Plots #
#########
percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}


#' Match summary plot
#'
#' @param results As produced by compare_alignments
#' @param ref Reference alignment
#' @param whether to display plot
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' PAC <- compare_alignments(ref,aln)
#' plot_match_summary(PAC)
plot_match_summary <- function(x,cys=FALSE,display=TRUE){
  
  identity       <- x$results[9,]
  proportion_cys <- (x$results[3,]/5)-0.2
  col = 1:ncol(x$results)
  plot_data = data.frame(Identity=identity,PropCys=proportion_cys,Column=col)
  
  p <- ggplot2::ggplot(plot_data,aes(x=Column)) + ggplot2::geom_line(aes(y=Identity,colour="% Identity"))
  if(cys) {
    p <- p + ggplot2::geom_line(aes(y=PropCys,colour="% Cysteines"))
    p <- p + ggplot2::geom_line(aes(y=0))
  }
  p <- p + ggplot2::theme(legend.title = element_blank())
  
  score <- x$score
  p <- p + ggplot2::geom_text(x=90,y=0.6,label=paste("Av =",percent(score)))
  
  if (display){
    print(p)
  }
  p
}


#' Category proportions plot
#'
#' @param results As produced by compare_alignments
#' @param whether to display plot
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' PAC <- compare_alignments(ref,aln)
#' plot_category_proportions(PAC)
plot_category_proportions <- function(x,stack=FALSE,display=TRUE){
  
  plot_data <- data.frame(Insertion=x$results[6,]/(1-x$results[5,]),
                          Deletion=x$results[7,]/(1-x$results[5,]),
                          Substitution=x$results[8,]/(1-x$results[5,]),
                          Column=1:ncol(x$results))
  md <- reshape2::melt(plot_data,id.vars='Column')
  colnames(md) <- c('Column','Change','Proportion')
  if (stack) {
    p <- ggplot2::ggplot(md,aes(x=Column,y=Proportion)) + ggplot2::geom_area(aes(fill=Change),position = 'stack')+ geom_line(aes(data=Change, ymax=1),position = 'stack')
  } 
  else {
    p <- ggplot2::ggplot(md,aes(x=Column,y=Proportion)) + ggplot2::geom_line(aes(color=Change))
  }
  
  if (display){
    print(p)
  }
  p
}

#' Alignment heatmap
#'
#' @param results As produced by compare_alignments
#' @param whether to display plot
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' PAC <- compare_alignments(ref,aln)
#' plot_alignment_heatmap(PAC)
plot_alignment_heatmap <- function(x,display=TRUE){
  
  hm_data <- t(t(x$means)/x$results[2,])
  md <- reshape2::melt(hm_data)
  colnames(md) <- c('Reference','Comparison','value')
  p <- ggplot2::ggplot(md) + ggplot2::geom_tile(aes(x=Reference,y=Comparison,fill=value)) + ggplot2::scale_fill_gradient("Agreement",low="white",high="black")
  p <- p  + ggplot2::theme(plot.background=element_rect(fill="white"),panel.background=element_rect(fill="white"))
  if (display){
    print(p)
  }
  p
}
#########









