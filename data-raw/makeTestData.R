###############
# Preparation #
###############

library("ade4")
library("seqinr")
percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

ref.fa <- "Alignment A.FA"
aln.fa <- "Alignment B.FA"

################
# FASTA import #
################

data.frame(read.fasta(ref.fa,set.attributes=FALSE)) -> ref   # convert to data frame of letters
data.frame(read.fasta(aln.fa,set.attributes=FALSE)) -> aln   # convert to data frame of letters

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
aln2    = matrix(nrow = dim(aln)[1], # number of aln columns
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
results = matrix(nrow = 10,           # matrix for results
                 ncol = dim(ref2)[1]) # number of ref2 columns
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

library(devtools)

devtools::use_data(ref2,overwrite =TRUE)
devtools::use_data(aln2,overwrite=TRUE)

# Calculating identity score
ci = 0
for(k in 1:dim(ref2)[1]){                         # for each (k) column of the ref2 alignment
  for(i in 1:dim(aln2)[1]){                       # for each (i) column of the aln2 alignment                               
    for(j in 1:dim(ref2)[2]){                     # for each (j) sequence
      ident[i,j] = identical(ref2[k,j],aln2[i,j])
      if ( ident[i,j] ){
        ci = ci + 1
        }# perform identity test
      means[i,k] = mean(ident[i,])                # all-against-all column identity (ref2 vs aln2)
    }
  }
  results[1,k] = which.max(means[,k])             # "ColumnMatch" which aln column had best match to each ref column
}
cat("CI ",ci)
devtools::use_data(results,overwrite=TRUE)
devtools::use_data(ident,overwrite=TRUE)

