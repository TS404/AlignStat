
# Alike Aligmnent Aligner

# Installation

```R
install.packages("devtools")
devtools::install_github("iracooke/AlikeAlignmentAligner")
```

# Example Usage

```R
data(ref)
data(aln)

results_list <- align_alignments(ref,aln)

alignment_heatmap(res_list$results,res_list$means,aln,ref)
```