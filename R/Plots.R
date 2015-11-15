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

  p <- ggplot(plot_data,aes(x=Column)) + geom_line(aes(y=Identity,colour="% Identity")) + geom_line(aes(y=-PropCys,colour="% Cysteines"))
  p <- p +  theme(legend.title = element_blank())

  score <- results[10,1]

  p <- p + geom_text(x=90,y=0.6,label=paste("Av =",percent(score)))
  if (display){
    print(p)
  }
  p
}

#' Category proportions plot
#'
#' @param results As produced by align_alignments
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' res_list <- align_alignments(ref,aln)
#' category_proportions_plot(res_list$results)
category_proportions_plot <- function(results){

  plot_data <- data.frame(Insertion=results[6,]/(1-results[5,]),
                         Deletion=results[7,]/(1-results[5,]),
                         Substitution=results[8,]/(1-results[5,]),
                         Column=1:ncol(results))
  md <- melt(plot_data,id.vars='Column')
  colnames(md) <- c('Column','Change','Proportion')
  p <- ggplot(md,aes(x=Column,y=Proportion)) + geom_line(aes(color=Change))
  
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
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' res_list <- align_alignments(ref,aln)
#' alignment_heatmap(res_list$results,res_list$means,aln,ref)
alignment_heatmap <- function(results,means,aln,ref){
  heatmap.2(t(t(means)/results[2,]),                            # gap-scaled data matrix
            trace        = "none",                              # no trace line
            key          = "FALSE",                             # no key
            dendrogram   = "none",                              # no dendrograms
            Rowv         = "FALSE",                             # no dendrograms
            Colv         = "FALSE",                             # no dendrograms
            density.info = "none",                              # no histogram
            labRow       = 1:dim(aln)[1]*c(rep(NA,9),TRUE),     # label every 10th row
            labCol       = 1:dim(ref)[1]*c(rep(NA,9),TRUE),     # label every 10th col
            xlab         = "Reference alignment",               # x axis label    
            ylab 	       = "Comparison alignment",              # y axis label
            col          = grey.colors(256, start=1, end=0))    # colours
}
#########









