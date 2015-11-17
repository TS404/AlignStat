Alike Aligmnent Aligner
=======================

>Thomas M A Shafee, Ira R Cooke, Marilyn A Anderson

>Department of Biochemistry, La Trobe Institute for Molecular Science, La Trobe University, Melbourne, Australia  
>College of Science, Health and Engineering, La Trobe University, Melbourne, Australia  
>Life Sciences Computation Centre, Victorian Life Sciences Computation Initiative, Melbourne, Australia

Title
-----
Comparison of alternative alignments of the same sequences to one another

Description
-----------
This package contains functions that compare two alternative multipe sequence alignemtns (MSAs) to determine whether they align homologous residues in the same columns as one another. It then classifies similarities and differences into conserved gaps, conserved sequence, insertion, deletion or substitution. Summarising these categories for each column yeilds information on which columns are agreed upon my both MSAs, and which differ. Several plotting functions easily visualise the comparison data for analysis.


Contains functions
------------------
```
align_alignments              Pairwise alignment of alignments
match_summary_plot            Summary plot of alignment similarities
category_proportions_plot     Detailed plot of alignment differences
alignment_heatmap             Heatmap of similarities between alignment columns
```

Installation
------------
```R
install.packages("devtools")
devtools::install_github("iracooke/AlikeAlignmentAligner")
```


align_alignments 
================
This function aligns two multiple sequence alignments (MSA) against one another. The alternative alignments must contain the same sequences in the same order. The function will classify any similarities and differences between the two MSAs. 

It is required as the first step any other package functions.

Description
-----------
Comparison of alternative multiple sequence alignments

###Usage
```
align_alignments (ref, aln,)
```

###Arguments
```
ref   The reference MSA (in fasta format)
aln   The MSA to compaer (in fasta format)
```

###Details
The `align_alignments` function first checks that the MSAs are alternative alignemnts of the same sequences. The function labels each character in the alignment with its occurence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns are the closest matches between the ref and aln MSAs. Each pairwise column comparison is stored as `something`

Each column is scores

###Example
```R
data(ref)
data(aln)
res_list <- align_alignments(ref,aln)
```







## Example Plots

```R
match_summary_plot       (res_list$results,ref)
category_proportions_plot (res_list$results)
alignment_heatmap         (res_list$results,res_list$means,aln,ref)
```
