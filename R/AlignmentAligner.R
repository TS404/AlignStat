#' Score the agreement between two alignments
#'
#' @param ref.fa Reference against which to score
#' @param aln.fa Alignment to score against reference
#' @return Results matrix with one column per sequence and the following rows
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
#' results <- align_alignments(ref,aln)
#'
align_alignments <- function(ref.fa,aln.fa){

  data.frame(read.fasta(ref.fa,set.attributes=FALSE)) -> ref
  data.frame(read.fasta(aln.fa,set.attributes=FALSE)) -> aln

  ###########################################
  # Replacing letters with letter+occurance #
  ###########################################

  # New matrix for numbered ref sequences
  ref.num = matrix(nrow = dim(ref)[1], # number of ref columns
                   ncol = dim(ref)[2]) # number sequences
  ref2    = matrix(nrow = dim(ref)[1], # number of ref columns
                   ncol = dim(ref)[2]) # number sequences

  # Number ref residues by their occurance
  for (s in 1:dim(ref)[2]){
    for (l in 1:dim(ref)[1]){
      ref.num[l,s] = sum(ref[1:l,s]==ref[l,s])
      paste(ref[,s],ref.num[,s])->ref2[,s]
    }  
  }

  # Remove extra space and de-number gaps
  gsub(x = ref2, pattern = " ",     replacement = "")  -> ref2
  gsub(x = ref2, pattern = "[-].*", replacement = "-") -> ref2

  # New matrix for numbered sequences
  aln.num = matrix(nrow = dim(aln)[1], # number of aln columns
                   ncol = dim(aln)[2]) # number sequences

  ##################

  # New matrix for numbered aln sequences
  aln.num = matrix(nrow = dim(aln)[1], # number of aln columns
                   ncol = dim(aln)[2]) # number sequences
  aln2 = matrix(nrow = dim(aln)[1], # number of aln columns
                  ncol = dim(aln)[2]) # number sequences

  # Number aln residues by their occurance
  for (s in 1:dim(aln)[2]){
    for (l in 1:dim(aln)[1]){
      aln.num[l,s] = sum(aln[1:l,s]==aln[l,s])
      paste(aln[,s],aln.num[,s])->aln2[,s]
    }  
  }

  # Remove extra space and de-number gaps
  gsub(x = aln2, pattern = " ",     replacement = "")  -> aln2
  gsub(x = aln2, pattern = "[-].*", replacement = "-") -> aln2

  # Replacing "-" with NA in the test alignment means
  # that gaps don't count towards column matching score
  aln[aln=="-"]  <-NA
  aln2[aln2=="-"]<-NA

  ##################################
  # Alignment identity calculation #
  ##################################

  ident = matrix(nrow = dim(aln2)[1],   # number of aln2 columns to test
                 ncol = dim(ref2)[2])   # number of sequences
  means = matrix(nrow = dim(aln2)[1],   # number of test alignment columns
                 ncol = dim(ref2)[1])   # number of ref2 columns

  # Making matrix for results data

  results = rcpp_align(ref2,aln2)
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

  ############################
  # Categorising differences # 
  ############################

  cat = matrix(nrow = dim(ref)[1], # number of ref columns
               ncol = dim(ref)[2]) # number sequences

  # Combining matrices
  for (x in 1:dim(ref)[1]){
    for (y in 1:dim(ref)[2]){
      paste(ref2[x,y],aln2[results[1,x],y])->cat[x,y]
    }
  }

  # Categorise mismatches
  cat2 <- cat
  gsub(x = cat2, pattern = "([a-z][0-9]+) \\1",       replacement = "M") -> cat2  # Match
  gsub(x = cat2, pattern = "[-] NA",                  replacement = "G") -> cat2  # Gap
  gsub(x = cat2, pattern = "[-] ([a-z][0-9]+)",       replacement = "I") -> cat2  # Insertion
  gsub(x = cat2, pattern = "([a-z][0-9]+) NA",        replacement = "D") -> cat2  # Deletion
  gsub(x = cat2, pattern = "[a-z][0-9]+ [a-z][0-9]+", replacement = "S") -> cat2  # Substitution

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
  results[10,1] <- sum(results[4,])/sum(1-results[5,])   # "Score"
  results[10,1] -> score
  results
}