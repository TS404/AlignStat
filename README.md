
# Alike Aligmnent Aligner

# Installation

install.packages("devtools")

devtools::install_github("iracooke/AlikeAlignmentAligner")

# Example Usage

data(ref)

data(aln)

res_list <- align_alignments(ref,aln)

alignment_heatmap(res_list$results,res_list$means,aln,ref)
