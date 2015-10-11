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
#' proportion_cys_plot(res_list$results,ref)
proportion_cys_plot <- function(results,ref){

  plot (results[9,],
        type = "l",
        ylim = c(-0.2, 1),
        xlab = "Column",
        ylab = "Identity")
  lines(-results[3,]/5, 
        type = "h")
  text (x = dim(ref)[1]+5,
        y = 0.9,
        label = paste("Av =",percent(score)),
        pos = 2)
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

  # Category proportions
  plot (results[6,]/(1-results[5,]), col="blue" ,type="l",
        ylim = c(-0, 1), xlab="Column", ylab="Category")         # Ins
  lines(results[7,]/(1-results[5,]), col="red"  ,type="l")       # Del
  lines(results[8,]/(1-results[5,]), col="green",type="l")       # Sub
  #lines(results[4,]/(1-results[5,]), col="black",type="l",lwd=2) # Match overlayed
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









