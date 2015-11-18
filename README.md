Alike Aligmnent Aligner
=======================

>Thomas M A Shafee, Ira R Cooke, Marilyn A Anderson

>Department of Biochemistry, La Trobe Institute for Molecular Science, La Trobe University, Melbourne, Australia  
>College of Science, Health and Engineering, La Trobe University, Melbourne, Australia  
>Life Sciences Computation Centre, Victorian Life Sciences Computation Initiative, Melbourne, Australia

Package title
-------------
Comparison of alternative multiple sequence alignments to one another

Description
-----------
This package contains functions that compare two alternative multipe sequence alignemtns (MSAs) to determine whether they align homologous residues in the same columns as one another. It then classifies similarities and differences into conserved gaps, conserved sequence, insertion, deletion or substitution. Summarising these categories for each column yeilds information on which columns are agreed upon my both MSAs, and which differ. Several plotting functions easily visualise the comparison data for analysis.


Contains functions
------------------
```R
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


----------------------------------------------------------------------------------------------
align_alignments 
================
This function aligns two multiple sequence alignments (MSA) against one another. The alternative alignments must contain the same sequences in the same order. The function will classify any similarities and differences between the two MSAs. 

It is required as the first step any other package functions.

Description
-----------
Comparison of alternative multiple sequence alignments

###Usage
```R
align_alignments (ref, aln)
```

###Arguments
```R
ref   The reference MSA (in fasta format)
aln   The MSA to compaer (in fasta format)
```

###Value
An object of class pairwise alignment comparison ('PAC') providing the optimal alignment of alignments and comparison of the differences between them.

The details of the output components are as follows:
```R
Results        A matrix with the following comparison statistics for each ith column
               of the reference alignment compared to its best match in the comparison
               alignemnt:

   Columnmatch   The column of the comparison alignment with the highest final match score
   Nongap        The proportion of characters that are not conserved gaps (all-Gapcon)
   Cys           The proportion of cysteines (relevant for cystieine rich proteins)
   Match         The proportion of characters that are identical between alignments
   Gapcon        The proportion of characters that are conserved gaps
   Insertion     The proportion of characters that are a gap in the reference, but
                 are a residue in the comparison alignment
   Deletion      The proportion of characters that are a residue in the reference,
                 but a gap in the comparison alignment
   Substitution  The proportion of characters that are one residue in the reference,
                 but but a non-homologous residue in the comparison alignment
   Finalmatch    The proportion of characters that match a a proportion of those
                 that are not conserved gaps (Match/Nongap)
   Score         The average Finalmatch score for all columns of the alignment 

Means      A matrix whose [i,k]th entry is the final match score between the ith
           column of the reference alignment and the kth column of the comparison
           alignment (Match/Nongap). Used to determine which columns are most similar
           for further analysis. Used to generate alignment heatmap
```

###Details
The `align_alignments` function first checks that the MSAs are alternative alignemnts of the same sequences. The function labels each character in the alignment with its occurence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns are the closest matches between the ref and aln MSAs. Each pairwise column comparison is stored as the `$means` value of the output. From this matrix, the comparison alignemtn column with the highest final match to each reference alignment column is used to calculate further statistics for the `$results` value of the output.

###Example
```R
data(ref)
data(aln)
res_list <- align_alignments(ref,aln)
```


---------------------------------------------------------------------------------------------
alignment_heatmap
=================
Description
-----------
Comparison of alternative multiple sequence alignments

Description
-----------
A heatmap of the column identities between two multiple sequence alignments

###Usage
```R
alignment_heatmap (results, means, aln, ref, display=TRUE)
```

###Arguments
```R
results     proportion of characters for each pairwise column comparison that are
            not conserved gaps (typically $results of results list generated by
            align_alignments)
means       proportion of identical characters for each pairwise column comparison
            (typically $means of results list generated by align_alignments)
ref         the reference alignment (to calculate consensus sequence)
aln         the comparison alignment (to calculate consensus sequence)
display     display this plot (default = TRUE)
```

###Details
The `alignment_heatmap` displays the similarity between each pairwise column comparison for the reference and comparison multiple sequence alignments. Colour density is determined by the proportion of identical character matches (M) between the columns, normalised to the number of characters that are not merely conserved gaps (not G). The (M) and (not G) values are generated by the `align_alignments` function as the result list items `$results` and `$means` respectively. This gives a representation of which columns are well agreed upon by the MSAs, and which columns are split by one MSA relative to the other.

###Example
```R
alignment_heatmap (res_list$results,res_list$means,aln,ref)
```


---------------------------------------------------------------------------------------------
match_summary_plot
=================

Description
-----------
A line plot summary of column similarity beween the two multiple sequence alignments 

###Usage
```R
match_summary_plot (results, cys=FALSE, display=TRUE)
```

###Arguments
```R
results    the pairwise alignment comparison results object 
           (typically the summary file generated by align_alignments)
cys        additionally show the cysteine abundance for each column (default = FALSE)
display    display this plot (default = TRUE)
```

###Details
This function generates a plot that sumarises the agreement between the two multiple sequence alignments. It plots the proportion of identical character matches (M) for each column as a proportion of the characteres that are not merely conserved gaps (not G). The overall average proportion of identical characters that are not conserved gaps is overlayed as a percentage. For alignments of cysteine-rich proteins, it can additionally plot the cysteine abundance for each column to indicate columns containing conserved cysteines (`cys=TRUE`).

###Example
```R
match_summary_plot (res_list$results,ref,cys=TRUE)
```


---------------------------------------------------------------------------------------------
category_proportions_plot
=================

Description
-----------
A line plot of the different causes of column dissimilarity beween the two multiple sequence alignments 

###Usage
```R
category_proportions_plot (results,stacked=FALSE,display=TRUE)
```

###Arguments
```R
results    the pairwise alignment comparison results object 
           (typically the summary file generated by align_alignments)
stacked    stacked area plot in stead of line plot (deault = FALSE)
display    display this plot (default = TRUE)
```

###Details
This function generates a detailed breakdown of the differences between the multiple sequence alignements. For all characters that are neither identical residues, nor conserved gaps, the relative proportions of insertions (I), deletions (D) and substitutions (S) is plotted.

###Example
```R
category_proportions_plot (res_list$results)
```

---------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------
Full example workflow
=====================

```R
# Alignment calculation
align_alignments          (ref, aln)
# Results visualisation
alignment_heatmap         (res_list$results,res_list$means,aln,ref)
match_summary_plot        (res_list$results,cys=TRUE)
category_proportions_plot (res_list$results,stacked=TRUE)
```
