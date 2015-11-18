#' Score the agreement between two alignments
#'
#' @param ref.fa Reference against which to score
#' @param aln.fa Alignment to score against reference
#' @return List containing results and means. Results is a matrix with one column per sequence and the following rows
#' \itemize{
#'  \item ColumnMatch 
#'  \item NonGap
#'  \item Cys
#'  \item Match
#'  \item Gap(con)
#'  \item Insertion
#'  \item Deletion
#'  \item Substitution
#'  \item FinalMatch
#'  \item Score
#' } 
#' 
#' @export
#' @examples
#' data(ref)
#' data(aln)
#' res_list <- compare_alignments(ref,aln)
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
                        "Gap(con)",     # 5
                        "Insertion",    # 6
                        "Deletion",     # 7
                        "Substitution", # 8
                        "FinalMatch",   # 9
                        "Score")        # 10
  
  cat2 <- res_list$cat
  
  # Write categories to results
  results[4,] <- t(rowMeans(cat2=="M")) # "Match"
  results[5,] <- t(rowMeans(cat2=="G")) # "Gap(con)"
  results[6,] <- t(rowMeans(cat2=="I")) # "Insertion"
  results[7,] <- t(rowMeans(cat2=="D")) # "Deletion"
  results[8,] <- t(rowMeans(cat2=="S")) # "Substitution"
  # Ref column gappiness
  results[2,] <- (1-t(rowMeans(ref=="-")))               # "NonGap" 
  # Ref cysteine occurance
  results[3,] <- (t(rowMeans(ref=="c")))                 # "Cys"
  # Correct Match scores to take gappiness into account
  results[9,] <- results[4,]/(1-results[5,])             # "FinalMatch"
  
  # Final mean!!!!!
  results[10,]  <- rep(NA,ncol(results))
  results[10,1] <- sum(results[4,])/sum(1-results[5,])   # "Score"
  # results[10,1] -> score
  
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
       alncon  = alncon)
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


