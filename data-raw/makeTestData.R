# Note: Run this script from package root not from inside raw_data


###############
# Preparation #
###############

library("ade4")
library("seqinr")
percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

ref.fa  <- "data-raw/Alignment A.fa"
com.aln <- "data-raw/Alignment B.clustal"

################
# FASTA import #
################

data.frame(read.fasta(ref.fa,set.attributes=FALSE)) -> ref   # convert to data frame of letters

temp <- seqinr::read.alignment(com.aln,format='clustal')
temp$nam <- do.call("rbind", lapply(strsplit(temp$nam," "),"[[", 1))
# reformat to data frame
output <- data.frame(strsplit(gsub("[\r\n]","",unlist(temp$seq)),split = ""))
colnames(output) <- temp$nam

output -> com   # convert to data frame of letters


reference_alignment <- ref
comparison_alignment <- com
devtools::use_data(reference_alignment,overwrite =TRUE)
devtools::use_data(comparison_alignment,overwrite=TRUE)


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
com.num = matrix(nrow = dim(com)[1], # number of com columns
                 ncol = dim(com)[2]) # number sequences

##################

# New matrix for numbered com sequences
com.num = matrix(nrow = dim(com)[1], # number of com columns
                 ncol = dim(com)[2]) # number sequences
com2    = matrix(nrow = dim(com)[1], # number of com columns
                 ncol = dim(com)[2]) # number sequences

# Number com residues by their occurance
for (s in 1:dim(com)[2]){
  for (l in 1:dim(com)[1]){
    com.num[l,s] = sum(com[1:l,s]==com[l,s])
    paste(com[,s],com.num[,s])->com2[,s]
  }  
}

# Remove extra space and de-number gaps
gsub(x = com2, pattern = " ",     replacement = "")  -> com2
gsub(x = com2, pattern = "[-].*", replacement = "-") -> com2

# Replacing "-" with NA in the test alignment means
# that gaps don't count towards column matching score
com[com=="-"]  <-NA
com2[com2=="-"]<-NA

##################################
# Alignment identity calculation #
##################################

ident = matrix(nrow = dim(com2)[1],   # number of com2 columns to test
               ncol = dim(ref2)[2])   # number of sequences
means = matrix(nrow = dim(com2)[1],   # number of test alignment columns
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

prepared_ref <- ref2
prepared_com <- com2

saveRDS(prepared_ref,file = "tests/testthat/prepared_ref.rda")
saveRDS(prepared_com,file = "tests/testthat/prepared_com.rda")
# devtools::use_data(prepared_ref,overwrite =TRUE)
# devtools::use_data(prepared_com,overwrite=TRUE)

# Calculating identity score
ci = 0
for(k in 1:dim(ref2)[1]){                         # for each (k) column of the ref2 alignment
  for(i in 1:dim(com2)[1]){                       # for each (i) column of the com2 alignment                               
    for(j in 1:dim(ref2)[2]){                     # for each (j) sequence
      ident[i,j] = identical(ref2[k,j],com2[i,j])
      if ( ident[i,j] ){
        ci = ci + 1
        }# perform identity test
      means[i,k] = mean(ident[i,])                # all-against-all column identity (ref2 vs com2)
    }
  }
  results[1,k] = which.max(means[,k])             # "ColumnMatch" which com column had best match to each ref column
}
cat("CI ",ci)
# devtools::use_data(means,overwrite=TRUE)
saveRDS(means,file = "tests/testthat/means.rda")

cat = matrix(nrow = dim(ref)[1], # number of ref columns
             ncol = dim(ref)[2]) # number sequences

# Combining matrices
for (x in 1:dim(ref)[1]){
  for (y in 1:dim(ref)[2]){
    paste(ref2[x,y],com2[results[1,x],y])->cat[x,y]
  }
}

# Categorise mismatches
cat2 <- cat
gsub(x = cat2, pattern = "([a-z][0-9]+) \\1",       replacement = "M") -> cat2  # Match
gsub(x = cat2, pattern = "[-] NA",                  replacement = "G") -> cat2  # Gap
gsub(x = cat2, pattern = "[-] ([a-z][0-9]+)",       replacement = "I") -> cat2  # Insertion
gsub(x = cat2, pattern = "([a-z][0-9]+) NA",        replacement = "D") -> cat2  # Deletion
gsub(x = cat2, pattern = "[a-z][0-9]+ [a-z][0-9]+", replacement = "S") -> cat2  # Substitution


categories <- cat2
# devtools::use_data(cat,overwrite=TRUE);
# devtools::use_data(categories,overwrite=TRUE);
saveRDS(categories,file = "tests/testthat/categories.rda")

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

# devtools::use_data(results,overwrite=TRUE)
saveRDS(results,file = "tests/testthat/results.rda")
