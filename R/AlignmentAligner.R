#' Compare alternative multiple sequence alignments
#'
#' @param ref   The reference MSA (in fasta format)
#' @param aln   The MSA to compare (in fasta format)
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
#' data(ref)
#' data(aln)
#' PAC <- compare_alignments(ref,aln)
#'
compare_alignments <- function(ref,aln){
  
  if (!is.data.frame(ref)){
    data.frame(seqinr::read.fasta(ref,set.attributes=FALSE)) -> ref
  } 
  if (!is.data.frame(aln)){
    data.frame(seqinr::read.fasta(aln,set.attributes=FALSE)) -> aln
  }
  
  if( !valid_alignments(ref,aln) ){
    stop("both alignments must contain the same sets of sequences in the same order")
  }
  
  
  ###########################################
  # Replacing letters with letter+occurance #
  ###########################################
  ref2 <- prepare_alignment_matrix(ref)
  aln2 <- prepare_alignment_matrix(aln)
  
  # Replacing "-" with NA in the test alignment means
  # that gaps don't count towards column matching score
  aln[aln=="-"]  <-NA
  aln2[aln2=="-"]<-NA
  
  ##################################
  # Alignment identity calculation #
  ##################################
  
  res_list = rcpp_align(ref2,aln2)
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
  alnlen <- nrow(aln)
  
  # Alignment consensus sequences
  refcon <- seqinr::consensus(t(ref))
  alncon <- seqinr::consensus(t(ref))
  
  # Create final object
  list(results = results,
       means   = res_list$means,
       cat     = cat,
       reflen  = reflen,
       alnlen  = alnlen,
       refcon  = refcon,
       alncon  = alncon,
       score   = score)
}


prepare_alignment_matrix <- function(alnmat){
  
  mat2 <- rcpp_prepare_alignment_matrix(as.matrix(alnmat))
  # Remove extra space and de-number gaps
  gsub(x = mat2, pattern = " ",     replacement = "")  -> mat2
  gsub(x = mat2, pattern = "[-].*", replacement = "-") -> mat2
  mat2
}


valid_alignments <- function(ref,aln){
  checks = sapply(1:ncol(ref),function(i){
    r=as.character(ref[,i])
    a=as.character(aln[,i])
    dg_ref = r[r!="-"]
    dg_aln = a[a!="-"]
    all(dg_aln==dg_ref)
  })
  all(checks)
}


