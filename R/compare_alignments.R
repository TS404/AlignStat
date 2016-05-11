#' Compare alternative multiple sequence alignments
#'
#' @param ref   The reference MSA (in fasta, clustal, msf, phylip or mase format)
#' @param com   The MSA to compare (in fasta, clustal, msf, phylip or mase format)
#' @param SP    Optionally also compute sum of pairs scores (default=FALSE)
#'
#' @return Generates an object of class "pairwise alignment comparison" (PAC), providing the optimal pairwise column alignment of two alternative MSAs of the same sequences, and summary statistics of the differences between them. The details of the PAC output components are as follows:
#' \itemize{
#'  \item {reference_P}          {The numbered character matrix of the reference alignment}
#'  \item {comparison_Q}         {The numbered character matrix of the comparison alignment}
#'  \item {results_R}            {The results summary matrix (containing column averages of match, gapcon, merge, split, shift)}
#'  \item {similarity_S}         {The similarity matrix between the reference and comparison alignment columns}
#'  \item {dissimilarity_D}      {The dissimilarity matrix between the reference and comparison (containing match, gapcon, merge, split, shift)}
#'  \item {dissimilarity_simple} {The dissimilarity matrix with categories stacked into a single 2D matrix}
#'  \item {columnmatch}          {The column of the comparison alignment with the highest final match score}
#'  \item {cys}                  {The proportion of cysteines (relevant for cysteine rich proteins)}
#'  \item {reflen}               {The length of the reference alignment}
#'  \item {comlen}               {The length of the comparison alignment}
#'  \item {refcon}               {The consensus sequence of the reference alignment}
#'  \item {comcon}               {The consensus sequence of the comparison alignment}
#'  \item {similarity_score}     {The overall similarity score}
#'  \item {sum_of_pairs}         {The sum of pairs score and related data (optional)}
#' } 
#' 
#' @export
#' @examples
#' data(reference_alignment)
#' data(comparison_alignment)
#' PAC <- compare_alignments(reference_alignment,comparison_alignment)
#'
#' @note The `compare_alignments` compares two alternative multiple sequence alignments (MSAs) of the same sequences. The alternative alignments must contain the same sequences. The function classifies similarities and differences between the two MSAs. It produces the "pairwise alignment comparison" object required as the first step any other package functions. The function converts the MSAs into matrices of sequence characters labelled by their occurrence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns have the highest similarty between the reference and comparison MSAs to generate a similarity matrix (excluding conserved gaps). From this matrix, the comparison alignment column with the similarity to each reference alignment column is used to calculate further statistics for dissimilarity matrix, summarised for each reference MSA column in the results matrix. Lastly, it calculates the overall similarity score between the two MSAs.
#'
compare_alignments <- function(ref,com,SP=FALSE){
  
  if (!is.data.frame(ref)){
    import_alignment(ref) -> ref
  }
  if (!is.data.frame(com)){
    import_alignment(com) -> com
  }
  # Check that all sequences are same present in both alignments even if different order
  if( !valid_alignments(ref,com)){
    stop("both alignments must contain the same sets of sequences, even if they are in a different order")
  }
  # Degap both alignments
  ref.degap <- degap_alignment(ref)
  com.degap <- degap_alignment(com)
  # Order full com alignment alphabetically
  com.alpha <- com[,order(com.degap)]
  # Reorder com.alpha by the order of ref.degap
  com <- com.alpha[,as.factor(ref.degap)]

  
  ###########################################
  # Replacing letters with letter+occurance #
  ###########################################
  ref2 <- prepare_alignment_matrix(ref)
  com2 <- prepare_alignment_matrix(com)
  
  if (!is.data.frame(ref)){
    names <- row.names(as.matrix(seqinr::read.fasta(ref)))
  } else {
    names <- colnames(ref)
  }
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
  catnames = c("Match","Gapcon","Merge","Split","Shift")
  
  row.names(results)<-c("ColumnMatch",  # 1
                        "NonGap",       # 2
                        "Cys",          # 3
                        "RawMatch",     # 4
                        "Gapcon",       # 5
                        "Merge",        # 6
                        "Split",        # 7
                        "Shift",        # 8
                        "Match")        # 9

  # Format P and Q matrices for output
  ref3 <- t(ref2)
  com3 <- t(com2)
  rownames(ref3) <- names
  rownames(com3) <- names
  # Restore "-"s to comparison alignment
  com3[is.na(com3)]<-"-"
  
  # Create dissimilarity (matrix D) from simplified dissimilarity (res_list$cat)
  dissimilarity_D <- array(dim      = c(ncol(ref), # rows
                                        nrow(ref), # columns
                                        5),        # stacks
                           dimnames = list(names,
                                           NULL,
                                           catnames))
                                           
  dissimilarity_D[,,1] <- 1*t(cat=="M") # "Match"
  dissimilarity_D[,,2] <- 1*t(cat=="g") # "Gapcon"
  dissimilarity_D[,,3] <- 1*t(cat=="m") # "Merge"
  dissimilarity_D[,,4] <- 1*t(cat=="s") # "Split"
  dissimilarity_D[,,5] <- 1*t(cat=="x") # "Shift"

  # Write category averages to results (R matrix)
  results_R           <- t(colMeans(dissimilarity_D))
  rownames(results_R) <- catnames

  # For each column of ref, which column of com is most similar
  columnmatch <- as.vector(res_list$results[1,]) 
 
  # Ref cysteine occurance
  cys <- as.vector((t(rowMeans(ref=="c"))))
 
  # Count alignment columns
  reflen <- nrow(ref)                              
  comlen <- nrow(com)
  
  # Alignment consensus sequences
  refcon <- seqinr::consensus(t(ref))
  comcon <- seqinr::consensus(t(com))

  # Final mean identity score
  similarity_score <- mean(cat=="M")/(1-mean(cat=="g"))
  
  ################
  # Sum of pairs #
  ################
  sum_of_pairs <- NULL
  
  if (SP==TRUE){
    
    P <- SPprep(ref3)
    ref.pairs     <- apply(t(P),1,list_pairs)
    ref.pairs.all <- unlist(ref.pairs)
    Q <- SPprep(com4)
    com.pairs     <- apply(t(Q),1,list_pairs)
    com.pairs.all <- unlist(com.pairs)
    
    SPSs <- NULL
    for(x in 1:length(com.pairs)){
      SPSs <- append(SPSs,SP(com.pairs[[x]],ref.pairs.all))
    }
    SPSs[is.nan(SPSs)] <- 0
    columnwise.SPS <- SPSs
    columnwise.CS  <- SPSs==1
    
    sum.of.pairs.score         <- SP(ref.pairs.all,com.pairs.all)
    reverse.sum.of.pairs.score <- PS(ref.pairs.all,com.pairs.all)
    column.score <- sum(columnwise.SPS==1)/comlen
    
    sum_of_pairs <- list(sum.of.pairs.score         = sum.of.pairs.score,
                         reverse.sum.of.pairs.score = reverse.sum.of.pairs.score,
                         columnwise.SPS             = columnwise.SPS,
                         column.score               = column.score,
                         columnwise.CS              = columnwise.CS)
  }


  # Create final object
  list(reference_P          = ref3,
       comparison_Q         = com4,
       results_R            = results_R,
       similarity_S         = means,
       dissimilarity_D      = dissimilarity_D,
       dissimilarity_simple = t(cat),
       columnmatch          = columnmatch,
       cys                  = cys,
       reflen               = reflen,
       comlen               = comlen,
       refcon               = refcon,
       comcon               = comcon,
       similarity_score     = similarity_score,
       sum_of_pairs         = sum_of_pairs)
}


import_alignment <- function(alignment,format=NULL){
  
  # default fmt
  fmt <- "fasta"
  
  # if clustal
  if( tools::file_ext(alignment)=="clustal"
     |tools::file_ext(alignment)=="CLUSTAL"
     |tools::file_ext(alignment)=="aln"
     |tools::file_ext(alignment)=="ALN"
     |tools::file_ext(alignment)=="clust"
     |tools::file_ext(alignment)=="clus"){
    fmt <- "clustal"
  }
  
  # if msf
  if( tools::file_ext(alignment)=="msf"
     |tools::file_ext(alignment)=="MSF"){
    fmt <- "msf"
  }
  
  # if mase
  if( tools::file_ext(alignment)=="mase"
     |tools::file_ext(alignment)=="MASE"){
    fmt <- "mase"
  }
  
  # if phylip
  if( tools::file_ext(alignment)=="phylip"
     |tools::file_ext(alignment)=="PHYLIP"
     |tools::file_ext(alignment)=="phy"
     |tools::file_ext(alignment)=="PHY"){
    fmt <- "phylip"
  }
  
  # format override
  if(!is.null(format)){
    fmt <- format
  }
  
  # import
  temp <- seqinr::read.alignment(alignment,format=fmt)
  # fix names
  temp$nam <- do.call("rbind", lapply(strsplit(temp$nam," "),"[[", 1))
  # reformat to data frame
  output <- data.frame(strsplit(unlist(temp$seq),split = ""))
  colnames(output) <- temp$nam
  
  output
}


prepare_alignment_matrix <- function(alignment){
  mat2 <- rcpp_prepare_alignment_matrix(as.matrix(alignment))
  # Remove extra space and de-number gaps
  gsub(x = mat2, pattern = " ",     replacement = "")  -> mat2
  gsub(x = mat2, pattern = "[-].*", replacement = "-") -> mat2
  mat2
}


degap_alignment <- function(final){
  # Remove gaps to convert alignment to list of strings
  gsub("-","",do.call("paste",c(data.frame(t(final)),sep="")))
}


valid_alignments <- function(ref,com){
  ref.degap <- degap_alignment(ref)
  com.degap <- degap_alignment(com)
  all(sort(ref.degap)==sort(com.degap))
}


# Fully unique identities fro all residues in MSA
SPprep <- function(x){
  matrix(paste(row.names(x),x,sep = "."),nrow = nrow(x))
}


# Full list of all pairs in an alignment column
list_pairs <- function(x){
  data <- x[-grep(pattern = "\\.[-]" , x)]
  tryCatch(do.call(paste,
                   as.data.frame(t(combn(data,2)),stringsAsFactors=FALSE)),
           error=function(e) NULL)
}


# Sum of pairs
SP <- function(reference,comparison){
  length(intersect(comparison,reference))/length(reference)
}


# Reverse sum of pairs
PS <- function(reference,comparison){
  length(intersect(comparison,reference))/length(comparison)
}
