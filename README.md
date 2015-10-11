
# Alike Aligmnent Aligner

# Installation

install.packages("devtools")
devtools::install_github("iracooke/AlikeAlignmentAligner")

# Usage

data(ref)
data(aln)

results <- align_alignments(ref,aln)