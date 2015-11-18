Alike Alignment Aligner
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
devtools::install_github("TS404/AlikeAlignmentAligner")
```


----------------------------------------------------------------------------------------------
compare_alignments 
================
This function aligns two multiple sequence alignments (MSA) against one another. The alternative alignments must contain the same sequences in the same order. The function will classify any similarities and differences between the two MSAs. 

It produces the "pairwise alignment comparison" object required as the first step any other package functions.

Description
-----------
Comparison of alternative multiple sequence alignments

###Usage
```R
compare_alignments (ref, aln)
```

###Arguments
```R
ref   The reference MSA (in fasta format)
aln   The MSA to compare (in fasta format)
```

###Value
Generates an object of class "pairwise alignment comparison" (PAC), providing the optimal alignment of alignments and comparison of the differences between them.

The details of the output components are as follows:
```R
results    A matrix with the following comparison statistics for each ith column of the
           reference alignment compared to its best match in the comparison alignment:

   Columnmatch   The column of the comparison alignment with the highest final match score
   Nongap        The proportion of characters that are not conserved gaps (all-Gapcon)
   Cys           The proportion of cysteines (relevant for cysteine rich proteins)
   Match         The proportion of characters that are identical between alignments
   Gapcon        The proportion of characters that are conserved gaps
   Insertion     The proportion of characters that are a gap in the reference, but are a
                 residue in the comparison alignment
   Deletion      The proportion of characters that are a residue in the reference, but a
                 gap in the comparison alignment
   Substitution  The proportion of characters that are one residue in the reference,
                 but a non-homologous residue in the comparison alignment
   Finalmatch    The proportion of characters that match a a proportion of those
                 that are not conserved gaps (Match/Nongap)
   Score         The average Finalmatch score for all columns of the alignment 

means      A matrix whose [i,k]th entry is the final match score between the ith
           column of the reference alignment and the kth column of the comparison
           alignment (Match/Nongap). Used to determine which columns are most similar
           for further analysis. Used to generate alignment heatmap

cat        A matrix whose [i,k]th entry is the match category of the ith sequence's
           kth residue for the reference alignment versus the comparison alignment
           (M=match, G=conserved gap, I=insertion, D=deletion, S=substitution)

reflen    number of columns in the reference alignment
alnlen    number of columns in the comparison alignment
refcon    consensus sequence of the reference alignment
alncon    consensus sequence of the comparison alignment
```

###Details
The `compare_alignments` function first checks that the MSAs are alternative alignments of the same sequences. The function labels each character in the alignment with its occurrence number in the sequence (e.g. to distinguish between the first and second cysteines of a sequence). It then compares the two MSAs to determine which columns are the closest matches between the ref and aln MSAs. Each pairwise column comparison is stored as the `$means` value of the output. From this matrix, the comparison alignment column with the highest final match to each reference alignment column is used to calculate further statistics for the `$results` value of the output.

###Example
```R
data(ref)
data(aln)
PAC <- compare_alignments(ref,aln)
```


---------------------------------------------------------------------------------------------
plot_alignment_heatmap
=================
Description
-----------
Comparison of alternative multiple sequence alignments

Description
-----------
A heatmap of the column identities between two multiple sequence alignments

###Usage
```R
plot_alignment_heatmap (results, means, aln, ref, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison"
           (typically the summary file generated by compare_alignments)
display    display this plot (default = TRUE)
```

###Details
The `plot_alignment_heatmap` displays the similarity between each pairwise column comparison for the reference and comparison multiple sequence alignments. Colour density is determined by the proportion of identical character matches (M) between the columns, normalised to the number of characters that are not merely conserved gaps (not G). The (M) and (not G) values are generated by the `compare_alignments` function as the result list items `$results` and `$means` respectively. This gives a representation of which columns are well agreed upon by the MSAs, and which columns are split by one MSA relative to the other.

###Example
```R
plot_alignment_heatmap (PAC,ref,aln)
```


---------------------------------------------------------------------------------------------
plot_match_summary
=================

Description
-----------
A line plot summary of column similarity between the two multiple sequence alignments 

###Usage
```R
plot_match_summary (results, cys=FALSE, display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison" 
           (typically the summary file generated by compare_alignments)
cys        additionally show the cysteine abundance for each column (default = FALSE)
display    display this plot (default = TRUE)
```

###Details
This function generates a plot that summarises the agreement between the two multiple sequence alignments. It plots the proportion of identical character matches (M) for each column as a proportion of the characters that are not merely conserved gaps (not G). The overall average proportion of identical characters that are not conserved gaps is overlaid as a percentage. For alignments of cysteine-rich proteins, it can additionally plot the cysteine abundance for each column to indicate columns containing conserved cysteines (`cys=TRUE`).

###Example
```R
plot_match_summary (PAC,cys=TRUE)
```


---------------------------------------------------------------------------------------------
plot_category_proportions
=================

Description
-----------
A line plot of the different causes of column dissimilarity between the two multiple sequence alignments 

###Usage
```R
plot_category_proportions (x,stack=FALSE,display=TRUE)
```

###Arguments
```R
x          an object of type "pairwise alignment comparison" 
           (typically the summary file generated by compare_alignments)
stack      stacked area plot in stead of line plot (deault = FALSE)
display    display this plot (default = TRUE)
```

###Details
This function generates a detailed breakdown of the differences between the multiple sequence alignments. For all characters that are neither identical residues, nor conserved gaps, the relative proportions of insertions (I), deletions (D) and substitutions (S) is plotted.

###Example
```R
plot_category_proportions (PAC)
```

---------------------------------------------------------------------------------------------
------------------------------------------------------------------------------------------------
Full example workflow
=====================

```R
# Load example data
data(ref)
data(aln)
# Alignment calculation
PAC <- compare_alignments (ref, aln)
# Results visualisation
plot_alignment_heatmap    (PAC)
plot_match_summary        (PAC, cys=TRUE)
plot_category_proportions (PAC, stack=TRUE)
```
