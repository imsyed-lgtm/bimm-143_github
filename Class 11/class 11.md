# Lab 11 quarto conversion
Iman Syed A18596789

\#Alpha Data Base (AFDB)

The EBI maintains the largest databaste of AlphaFold structure
prediction models at: <http://alphafold.ebi.ac.uk> From last class
(before Halloween) we saw that the PDB had 244,290 (Oct 2025)

the total number of protien sequences in UniProtKB is 199,579,901

> Key Point: This is a tiny fraction of sequence space that has
> structural coverarge (0.12%)

``` r
(244290/199579901) * 100
```

    [1] 0.1224021

AFDB is attempting to address this gap…

There are two “quality scores” from AlphaFold- one for the residues(i.e
for each amino acid) called the **pLDDT** score. The other **PAE** score
that measures the confidence in the relative position of two residues
(i.e a score for every pair of residues)

\##Generating your own structure predictions

Figure of 5 generated HIV-PR models:

![](HIV.png)

``` r
#con <- consenus (aln,cutoff =0.9)
#con$seq
```

## Custom analysis of resulting models

``` r
# Change this for YOUR results dir name
results_dir <- "HIVPrDimertest2_23119" 

# File names for all PDB models
pdb_files <- list.files(path=results_dir,
                        pattern="*.pdb",
                        full.names = TRUE)

# Print our PDB file names
basename(pdb_files)
```

    character(0)

``` r
library(bio3d)

# Read all data from Models 
#  and superpose/fit coords


pdb_files <- list.files(pattern = ".pdb")
```
