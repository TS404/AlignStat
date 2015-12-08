#' AlignStat: A tool for the statistical comparison of alternative multiple sequence alignments
#' 
#' This package contains functions that compare two alternative multiple sequence 
#' alignments (MSAs) to determine whether they align homologous residues in the 
#' same columns as one another. It then classifies similarities and differences into 
#' conserved gaps, conserved sequence, insertion, deletion or substitution. 
#' Summarising these categories for each column yields information on which columns 
#' are agreed upon my both MSAs, and which differ. Several plotting functions easily 
#' visualise the comparison data for analysis.
#' 
#' @section Computing statistics:
#' 
#' Use \code{\link{compare_alignments}} to calculate stats on a pair of alignments
#'
#' @section Plotting functions:

#' Use \code{\link{plot_alignment_heatmap}} to view a heatmap of agreement
#' Use \code{\link{plot_match_summary}} to visualise agreement and cysteine content
#' Use \code{\link{plot_category_proportions}} to see various metrics of agreement
#' 
#' 
#' @docType package
#' @name AlignStat
#' @useDynLib AlignStat
#' @importFrom Rcpp sourceCpp
NULL