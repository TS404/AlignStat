
# Alike Aligmnent Aligner

## Installation

```R
install.packages("devtools")
devtools::install_github("iracooke/AlikeAlignmentAligner")
```

## Example Usage

```R
data(ref)
data(aln)

res_list <- align_alignments(ref,aln)
```

## Example Plots

```R
proportion_cys_plot       (res_list$results,ref)
category_proportions_plot (res_list$results)
alignment_heatmap         (res_list$results,res_list$means,aln,ref)
```
