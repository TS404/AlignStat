AlignStat: A tool for the statistical comparison of alternative multiple sequence alignments
=======================

>Thomas M A Shafee, Ira R Cooke, Marilyn A Anderson

>Department of Biochemistry, La Trobe Institute for Molecular Science, La Trobe University, Melbourne, Australia  
>College of Science, Health and Engineering, La Trobe University, Melbourne, Australia  
>Life Sciences Computation Centre, Victorian Life Sciences Computation Initiative, Melbourne, Australia

Package title
-------------
Comparison of alternative multiple sequence alignments

Description
-----------
This package contains functions that compare two alternative multiple sequence alignments (MSAs) to determine whether they align homologous residues in the same columns as one another. It then classifies similarities and differences into conserved gaps, conserved sequence, insertion, deletion or substitution. Summarising these categories for each column yields information on which columns are agreed upon my both MSAs, and which differ. Several plotting functions easily visualise the comparison data for analysis.


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
Generates an object of class "pairwise alignment comparison" (PAC), providing the optimal alignment of alignments and comparison of the differences between them. The details of the PAC output components are as follows:
```R
results    A matrix with the following comparison statistics for each ith column of the
           reference alignment compared to its best match in the comparison alignment:

   Columnmatch   The column of the comparison alignment with the highest final match score
   Cys           The proportion of cysteines (relevant for cysteine rich proteins)
   Match         The proportion of characters that are identical between alignments
   Gapcon        The proportion of characters that are conserved gaps
   Insertion     The proportion of characters that are a gap in the reference, but are a
                 residue in the comparison alignment
   Deletion      The proportion of characters that are a residue in the reference, but a
                 gap in the comparison alignment
   Substitution  The proportion of characters that are one residue in the reference,
                 but a non-homologous residue in the comparison alignment
   Finalmatch    The proportion of characters that match as a proportion of those
                 that are not conserved gaps (Match/not_Gapcon)

means      A matrix whose [i,k]th entry is the final match score between the ith
           column of the reference alignment and the kth column of the comparison
           alignment (Match/not_Gapcon). Used to determine which columns are most similar
           for further analysis. Used to generate alignment heatmap

cat        A matrix whose [i,k]th entry is the match category kth residue of the
           ith sequence for the reference alignment versus the comparison alignment
           (M=match, G=gapcon, I=insertion, D=deletion, S=substitution)

reflen    The number of columns in the reference alignment
comlen    The number of columns in the comparison alignment
refcon    Consensus sequence of the reference alignment
comcon    Consensus sequence of the comparison alignment

score     The average Finalmatch score for all columns of the alignment
```

###Details
The `compare_alignments` function first checks that the MSAs are alternative alignments of the same sequences. The function labels each character in the alignment with its occurrence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns are the closest matches between the ref and com MSAs. Each pairwise column comparison is stored as the `$means` value of the output. From this matrix, the comparison alignment column with the highest final match to each reference alignment column is used to calculate further statistics for the `$results` value of the output. The overall Finalmatch score for the whole comparison is output as the `$score` value.

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
plot_similarity_heatmap (x, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison"
           (typically the summary file generated by compare_alignments)
display    display this plot (default = TRUE)
```

###Details
This function displays the similarity between each pairwise column comparison for the reference and comparison MSAs. Colour density is determined by the proportion of identical character matches between the columns, normalised to the number of characters that are not merely conserved gaps. This gives a representation of which columns are well agreed upon by the MSAs, and which columns are split by one MSA relative to the other.

###Example
```R
plot_similarity_heatmap (PAC)
```


---------------------------------------------------------------------------------------------
plot_dissimilarity_matrix
=================
Plot a heatmap of the dissimilarity matrix of two multiple sequence alignments

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
This function displays the dissimilarity categories for all characters in the reference alignment. This gives a representation of which columns are well agreed upon by the MSAs, and which sequence regions of the reference alignment are split, merged, or shifted.

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
plot_similarity_summary (results, cys=FALSE, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison" 
           (typically the summary file generated by compare_alignments)
cys        additionally show the cysteine abundance for each column (default = FALSE)
display    display this plot (default = TRUE)
```

###Details
This function generates a plot that summarises the similarity between the two multiple sequence alignments for each column of the reference alignment. For each column, it plots the proportion of identical character matches as a proportion of the characters that are not merely conserved gaps. The overall average proportion of identical characters that are not conserved gaps is overlaid as a percentage. For alignments of cysteine-rich proteins, the cysteine abundance for each column may also be plotted to indicate columns containing conserved cysteines (`cys=TRUE`).

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
plot_category_proportions (x,stack=TRUE,display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison" 
           (typically the summary file generated by compare_alignments)
stack      stacked area plot in stead of line plot (default = TRUE)
display    display this plot (default = TRUE)
```

###Details
This function generates a detailed breakdown of the differences between the multiple sequence alignments for each column of the reference alignment. For each column, the relative proportions of merges, splits and shifts is plotted as a proportion of characters that are not merely conserved gaps.

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
