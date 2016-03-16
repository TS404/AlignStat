AlignStat: A tool for the statistical comparison of alternative multiple sequence alignments
=======================

>Thomas M A Shafee, Ira R Cooke

>Department of Biochemistry, La Trobe Institute for Molecular Science, La Trobe University, Melbourne, Australia  
>College of Science, Health and Engineering, La Trobe University, Melbourne, Australia  
>Life Sciences Computation Centre, Victorian Life Sciences Computation Initiative, Melbourne, Australia

Package title
-------------
Comparison of alternative multiple sequence alignments

Description
-----------
This package contains functions that compare two alternative multiple sequence alignments (MSAs) to determine how well they align homologous residues in the same columns as one another. It classifies similarities and differences into conserved sequence, conserved gaps, splits, merges and shifts. Summarising these categories for each column yields information on which columns are agreed upon by both MSAs, and which differ. Output graphs visualise the comparison data for analysis.


Contains functions
------------------
```R
compare_alignments            Pairwise alignment of alignments
plot_match_summary            Summary plot of alignment similarities
plot_category_proportions     Detailed plot of alignment differences
plot_alignment_heatmap        Heatmap of similarities between alignment columns
```

Installation
------------
From CRAN
```R
install.packages("AlignStat")
```
From GitHub
```R
install.packages("devtools")
devtools::install_github("TS404/AlignStat")
library("AlignStat")
```


----------------------------------------------------------------------------------------------
compare_alignments 
================
Compare alternative multiple sequence alignments

Description
-----------
This function aligns two multiple sequence alignments (MSA) against one another. The alternative alignments must contain the same sequences in the same order. The function will classify any similarities and differences between the two MSAs. 

It produces the "pairwise alignment comparison" object required as the first step any of the other package functions.


###Usage
```R
compare_alignments (ref, com)
```

###Arguments
```R
ref   The reference MSA (in fasta format)
com   The MSA to compare (in fasta format)
```

###Value
Generates an object of class "pairwise alignment comparison" (PAC), providing the optimal pairwise column alignment of two alternative MSAs of the same sequences, and summary statistics of the differences between them. The details of the PAC output components are as follows:
```R
reference_P           A numbered character matrix of the reference alignment
comparison_Q          A numbered character matrix of the comparison alignment
results_R             A matrix whos [i,j]th entry is the ith match category average of the
                      jth column of the reference alignment versus the comparison alignment
                      (i1=match, i2=conserved gap, i3=merge, i4=split, i5=shift) Used to
                      generate the similarity summary and dissimilarity summary plots.
similarity_S          A similarity matrix whose [i,j]th entry is the similarity score between
                      the ith column of the reference alignment and the jth column of the
                      comparison alignment. Used to determine which columns are most similar
                      for further analysis. Used to generate the similarity heatmap plot.
dissimilarity_D       A dissimilarity matrix whose [i,j,k]th entry is the kth match category
                      of the jth residue of the ith sequence for the reference alignment
                      versus the comparison alignment (k1=match, k2=conserved gap, k3=merge,
                      k4=split, k5=shift).
dissimilarity_simple  A matrix whose [i,j]th entry is the dissimilarity category of the jth
                      residue of the ith sequence for the reference alignment versus the
                      comparison alignment (M=match, g=conserved gap, m=merge, s=split, x=shift).
                      Generated from the dissimilarity matrix with categories stacked into a
                      single 2D matrix. Used to the dissimilarity matrix plot.
columnmatch           The column of the comparison alignment with the highest final match score
cys                   The proportion of cysteines (relevant for cysteine rich proteins)
reflen                The number of columns in the reference alignment
comlen                The number of columns in the comparison alignment
refcon                The consensus sequence of the reference alignment
comcon                The consensus sequence of the comparison alignment
score                 The overall similarity score between the reference and comparison alignments
```

###Details

The `compare_alignments` compares two alternative multiple sequence alignments (MSAs) of the same sequences. The alternative alignments must contain the same sequences in the same order. The function classifies similarities and differences between the two MSAs. It produces the "pairwise alignment comparison" object required as the first step any other package functions.

The function converts the MSAs into matrices of sequence characters labelled by their occurrence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns have the highest similarty between the reference and comparison MSAs to generate a similarity matrix (excluding conserved gaps). From this matrix, the comparison alignment column with the similarity to each reference alignment column is used to calculate further statistics for dissimilarity matrix, summarised for each reference MSA column in the results matrix. Lastly, it calculates the overall similarity score between the two MSAs.

###Example
```R
data("reference_alignment")
data("comparison_alignment")
PAC <- compare_alignments(reference_alignment,comparison_alignment)
```


---------------------------------------------------------------------------------------------
plot_similarity_heatmap
=================

A heatmap plot of the column identities between two multiple sequence alignments

###Usage
```R
plot_similarity_heatmap (x, scale=TRUE, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison"
           (typically the summary file generated by compare_alignments)
scale      scale data to proportion of characters that are not conserved gaps (default = TRUE)
display    display this plot (default = TRUE)
```

###Details
The `plot_similarity_heatmap` function displays the similarity between each pairwise column comparison for the reference and comparison MSAs. Colour density is determined by the proportion of identical character matches between the columns, normalised to the number of characters that are not merely conserved gaps. This gives a representation of which columns are well agreed upon by the MSAs, and which columns are split by one MSA relative to the other.

###Example
```R
plot_similarity_heatmap (PAC)
```


---------------------------------------------------------------------------------------------
plot_dissimilarity_matrix
=================
A heatmap plot of the dissimilarity matrix of two multiple sequence alignments

###Usage
```R
plot_dissimilarity_matrix (x, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison" 
           (typically the summary file generated by compare_alignments)
display    display this plot (default = TRUE)
```

###Details
The `plot_dissimilarity_matrix` function displays the dissimilarity categories for all characters in the reference alignment. This gives a representation of which columns are well agreed upon by the MSAs, and which sequence regions of the reference alignment are split, merged, or shifted.

###Example
```R
data(reference_alignment)
data(comparison_alignment)
PAC <- compare_alignments(reference_alignment,comparison_alignment)
plot_dissimilarity_matrix(PAC)
```


---------------------------------------------------------------------------------------------
plot_similarity_summary
=================

A line plot summary of column similarity between two multiple sequence alignments 

###Usage
```R
plot_similarity_summary (results, scale=TRUE, cys=FALSE, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison" 
           (typically the summary file generated by compare_alignments)
scale      scale data to proportion of characters that are not conserved gaps (default = TRUE)
cys        additionally show the cysteine abundance for each column (default = FALSE)
display    display this plot (default = TRUE)
```

###Details
The `plot_similarity_summary` function generates a plot that summarises the similarity between the two multiple sequence alignments for each column of the reference alignment. For each column, it plots the proportion of identical character matches as a proportion of the characters that are not merely conserved gaps. The overall average proportion of identical characters that are not conserved gaps is overlaid as a percentage. For alignments of cysteine-rich proteins, the cysteine abundance for each column may also be plotted to indicate columns containing conserved cysteines (`cys=TRUE`).

###Example
```R
plot_similarity_summary (PAC, cys=TRUE)
```


---------------------------------------------------------------------------------------------
plot_dissimilarity_summary
=================

An area plot summary of the different causes of column dissimilarity between two multiple sequence alignments

###Usage
```R
plot_category_proportions (x, scale=TRUE, stack=TRUE, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison" 
           (typically the summary file generated by compare_alignments)
scale      scale data to proportion of characters that are not conserved gaps (default = TRUE)
stack      stacked area plot in stead of line plot (default = TRUE)
display    display this plot (default = TRUE)
```

###Details
The `plot_dissimilarity_summary` function generates a detailed breakdown of the differences between the multiple sequence alignments for each column of the reference alignment. For each column, the relative proportions of merges, splits and shifts is plotted as a proportion of characters that are not merely conserved gaps.

###Example
```R
plot_dissimilarity_summary (PAC)
```

------------------------------------------------------------------------------------------------
Full example workflow
=====================

```R
# Example data loading
data("reference_alignment")
data("comparison_alignment")
# Alignment comparison calculation
PAC <- compare_alignments (reference_alignment,comparison_alignment)
# Results visualisation
plot_similarity_heatmap    (PAC)
plot_dissimilarity_matrix  (PAC)
plot_similarity_summary    (PAC, cys=TRUE)
plot_dissimilarity_summary (PAC, stack=TRUE)
```
