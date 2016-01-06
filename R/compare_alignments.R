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
#'  \item {reference_P}          {The numbered character matrix of the reference alignment}
#'  \item {comparison_Q}         {The numbered character matrix of the comparison alignment}
#'  \item {results_R}            {The results summary matrix (containing column averages of match, merge, split, shift, gapcon)}
#'  \item {similarity_S}         {The similarity matrix between the reference and comparison alignment columns}
#'  \item {dissimilarity_D}      {The dissimilarity matrix between the reference and comparison (containing match, merge, split, shift, gapcon)}
#'  \item {dissimilarity_simple} {The dissimilarity matrix with categories stacked into a single 2D matrix}
#'  \item {columnmatch}          {The column of the comparison alignment with the highest final match score}
#'  \item {cys}                  {The proportion of cysteines (relevant for cysteine rich proteins)}
#'  \item {reflen}               {The length of the reference alignment}
#'  \item {comlen}               {The length of the comparison alignment}
#'  \item {refcon}               {The consensus sequence of the reference alignment}
#'  \item {comcon}               {The consensus sequence of the comparison alignment}
#'  \item {score}                {The overall similarity score}
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
  
  # res_list contains $results, $cat, $means 
  res_list = rcpp_align(ref2,com2)
  results  = res_list$results
  cat      = res_list$cat
  means    = res_list$means
  
  row.names(results)<-c("ColumnMatch",  # 1
                        "NonGap",       # 2
                        "Cys",          # 3
                        "RawMatch",     # 4
                        "Gapcon",       # 5
                        "Merge",        # 6
                        "Split",        # 7
                        "Shift",        # 8
                        "Match")        # 9

  # Create dissimilarity (matrix D) from simplified dissimilarity (res_list$cat)
  dissimilarity_D <- array(, dim=c(ncol(ref), # rows
                                   nrow(ref), # columns
                                   5))        # stacks
  dissimilarity_D[,,1] <- 1*(cat=="M") # "Match"
  dissimilarity_D[,,2] <- 1*(cat=="m") # "Merge"
  dissimilarity_D[,,3] <- 1*(cat=="s") # "Split"
  dissimilarity_D[,,4] <- 1*(cat=="x") # "Shift"
  dissimilarity_D[,,5] <- 1*(cat=="g") # "Gapcon"

  # Write category averages to results (R matrix)
  results_R <- array(, dim=c(5,          # rows
                             nrow(ref))) # columns
  results_R[1,] <- t(rowMeans(cat=="M")) # "Match"
  results_R[2,] <- t(rowMeans(cat=="m")) # "Merge"
  results_R[3,] <- t(rowMeans(cat=="s")) # "Split"
  results_R[4,] <- t(rowMeans(cat=="x")) # "Shift"
  results_R[5,] <- t(rowMeans(cat=="g")) # "Gapcon"

  # For each column of ref, which column of com is most similar
  columnmatch <- res_list$results[1,] 
 
  # Ref cysteine occurance
  cys <- (t(rowMeans(ref=="c")))        
 
  # Count alignment columns
  reflen <- nrow(ref)                              
  comlen <- nrow(com)
  
  # Alignment consensus sequences
  refcon <- seqinr::consensus(t(ref))
  comcon <- seqinr::consensus(t(com))

  # Final mean score
  score <- mean(cat=="M")/(1-mean(cat=="g"))
  
  # Create final object
  list(reference_P          = t(ref2),
       comparison_Q         = t(com2),
       results_R            = results_R,
       similiarity_S        = means,
       dissimilarity_D      = dissimilarity_D,
       dissimilarity_simple = cat,
       columnmatch          = columnmatch,
       cys                  = cys,
       reflen               = reflen,
       comlen               = comlen,
       refcon               = refcon,
       comcon               = comcon,
       score                = score)
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
