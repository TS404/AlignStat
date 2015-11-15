#########
# Plots #
#########
percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}


#' Plot match proportion + Cys
#'
#' @param results As produced by align_alignments
#' @param ref Reference alignment
#' @param whether to display plot
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' res_list <- align_alignments(ref,aln)
#' match_summary_plot(res_list$results,ref)
match_summary_plot <- function(results,ref,display=TRUE){

  identity <- results[9,]
  proportion_cys <- results[3,]/5
  col = 1:ncol(results)
  plot_data = data.frame(Identity=identity,PropCys=proportion_cys,Column=col)

  p <- ggplot2::ggplot(plot_data,aes(x=Column)) + ggplot2::geom_line(aes(y=Identity,colour="% Identity")) + ggplot2::geom_line(aes(y=-PropCys,colour="% Cysteines"))
  p <- p +  ggplot2::theme(legend.title = element_blank())

  score <- results[10,1]

  p <- p + ggplot2::geom_text(x=90,y=0.6,label=paste("Av =",percent(score)))
  if (display){
    print(p)
  }
  p
}

#' Category proportions plot
#'
#' @param results As produced by align_alignments
#' @param whether to display plot
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' res_list <- align_alignments(ref,aln)
#' category_proportions_plot(res_list$results)
category_proportions_plot <- function(results,display=TRUE){

  plot_data <- data.frame(Insertion=results[6,]/(1-results[5,]),
                         Deletion=results[7,]/(1-results[5,]),
                         Substitution=results[8,]/(1-results[5,]),
                         Column=1:ncol(results))
  md <- reshape2::melt(plot_data,id.vars='Column')
  colnames(md) <- c('Column','Change','Proportion')
  p <- ggplot2::ggplot(md,aes(x=Column,y=Proportion)) + ggplot2::geom_line(aes(color=Change))
  
  if (display){
    print(p)
  }
  p
}

#' Alignment heatmap
#'
#' @param results As produced by align_alignments
#' @param aln Alignment to align to ref
#' @param ref Reference alignment
#' @param whether to display plot
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' res_list <- align_alignments(ref,aln)
#' alignment_heatmap(res_list$results,res_list$means,aln,ref)
alignment_heatmap <- function(results,means,aln,ref,display=TRUE){
  
  hm_data <- t(t(means)/results[2,])
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









