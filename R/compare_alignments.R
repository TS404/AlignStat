#' Compare alternative multiple sequence alignments
#'
#' @note This function aligns two multiple sequence alignments (MSA) against one another. The alternative alignments must contain the same sequences in the same order. The function will classify any similarities and differences between the two MSAs. It produces the "pairwise alignment comparison" object required as the first step any other package functions.
#' @note This The `compare_alignments` function first checks that the MSAs are alternative alignments of the same sequences. The function labels each character in the alignment with its occurrence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns are the closest matches between the ref and com MSAs. Each pairwise column comparison is stored as the `$means` value of the output. From this matrix, the comparison alignment column with the highest final match to each reference alignment column is used to calculate further statistics for the `$results` value of the output. The overall Finalmatch score for the whole comparison is output as the `$score` value.
#'
#' @param ref   The reference MSA (in fasta format)
#' @param com   The MSA to compare (in fasta format)
#'
#' @return Generates an object of class "pairwise alignment comparison" (PAC), providing the optimal alignment of alignments and comparison of the differences between them. The details of the PAC output components are as follows:
#' \itemize{
#'  \item {Columnmatch}   {The column of the comparison alignment with the highest final match score}
#'  \item {Nongap}        {The proportion of characters that are not gaps}
#'  \item {Cys}           {The proportion of cysteines (relevant for cysteine rich proteins)}
#'  \item {Match}         {The proportion of characters that are identical between alignments}
#'  \item {Gapcon}        {The proportion of characters that are conserved gaps}
#'  \item {Insertion}     {The proportion of characters that are a gap in the reference, but are a residue in the comparison alignment}
#'  \item {Deletion}      {The proportion of characters that are a residue in the reference, but a gap in the comparison alignment}
#'  \item {Substitution}  {The proportion of characters that are one residue in the reference, but a non-homologous residue in the comparison alignment}
#'  \item {Finalmatch}    {The proportion of characters that match as a proportion of those that are not conserved gaps (Match/not_Gapcon)}
#' } 
#' 
#' @export
#' @examples
#' data(reference_alignment)
#' data(comparison_alignment)
#' PAC <- compare_alignments(reference_alignment,comparison_alignment)
#'
compare_alignments <- function(ref,com){
  
  if (!is.data.frame(ref)){
    data.frame(seqinr::read.fasta(ref,set.attributes=FALSE)) -> ref
  } 
  if (!is.data.frame(com)){
    data.frame(seqinr::read.fasta(com,set.attributes=FALSE)) -> com
  }
  
  if( !valid_alignments(ref,com) ){
    stop("both alignments must contain the same sets of sequences in the same order")
  }
  
  
  ###########################################
  # Replacing letters with letter+occurance #
  ###########################################
  ref2 <- prepare_alignment_matrix(ref)
  com2 <- prepare_alignment_matrix(com)
  
  # Replacing "-" with NA in the test alignment means
  # that gaps don't count towards column matching score
  com[com=="-"]  <-NA
  com2[com2=="-"]<-NA
  
  ##################################
  # Alignment identity calculation #
  ##################################
  
  res_list = rcpp_align(ref2,com2)
  results  = res_list$results
  row.names(results)<-c("ColumnMatch",  # 1
                        "NonGap",       # 2
                        "Cys",          # 3
                        "Match",        # 4
                        "Gapcon",       # 5
                        "Insertion",    # 6
                        "Deletion",     # 7
                        "Substitution", # 8
                        "FinalMatch")   # 9
  
  cat2 <- res_list$cat
  
  # Write categories to results
  results[4,] <- t(rowMeans(cat2=="M")) # "Match"
  results[5,] <- t(rowMeans(cat2=="G")) # "Gapcon"
  results[6,] <- t(rowMeans(cat2=="I")) # "Insertion"
  results[7,] <- t(rowMeans(cat2=="D")) # "Deletion"
  results[8,] <- t(rowMeans(cat2=="S")) # "Substitution"
  # Ref column gappiness
  results[2,] <- (1-t(rowMeans(ref=="-")))               # "NonGap" 
  # Ref cysteine occurance
  results[3,] <- (t(rowMeans(ref=="c")))                 # "Cys"
  # Correct Match scores to take gappiness into account
  results[9,] <- results[4,]/(1-results[5,])             # "FinalMatch"
  
  # Final mean score
  score <- sum(results[4,])/sum(1-results[5,])           # "Score"
  
  # Count alignment columns
  reflen <- nrow(ref)
  comlen <- nrow(com)
  
  # Alignment consensus sequences
  refcon <- seqinr::consensus(t(ref))
  comcon <- seqinr::consensus(t(ref))
  
  # Create final object
  list(results = results,
       means   = res_list$means,
       cat     = cat,
       reflen  = reflen,
       comlen  = comlen,
       refcon  = refcon,
       comcon  = comcon,
       score   = score)
}


prepare_alignment_matrix <- function(commat){
  
  mat2 <- rcpp_prepare_alignment_matrix(as.matrix(commat))
  # Remove extra space and de-number gaps
  gsub(x = mat2, pattern = " ",     replacement = "")  -> mat2
  gsub(x = mat2, pattern = "[-].*", replacement = "-") -> mat2
  mat2
}


valid_alignments <- function(ref,com){
  checks = sapply(1:ncol(ref),function(i){
    r=as.character(ref[,i])
    a=as.character(com[,i])
    dg_ref = r[r!="-"]
    dg_com = a[a!="-"]
    all(dg_com==dg_ref)
  })
  all(checks)
}

