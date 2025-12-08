# Class 12 quarto


## 

##Background today we will analyze some RNASeq from Himes et.al on the
effects of a common steriod(dexamethasone) on airway smooth mucsle
cells.

Our starting point is the “counts” data and meta data that contain the
count values for each gene in their different experiements.

##Data import

``` r
# Complete the missing code
counts <- read.csv("airway_scaledcounts.csv", row.names=1)
metadata <-  read.csv("airway_metadata.csv")
```

``` r
head(counts)
```

                    SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516
    ENSG00000000003        723        486        904        445       1170
    ENSG00000000005          0          0          0          0          0
    ENSG00000000419        467        523        616        371        582
    ENSG00000000457        347        258        364        237        318
    ENSG00000000460         96         81         73         66        118
    ENSG00000000938          0          0          1          0          2
                    SRR1039517 SRR1039520 SRR1039521
    ENSG00000000003       1097        806        604
    ENSG00000000005          0          0          0
    ENSG00000000419        781        417        509
    ENSG00000000457        447        330        324
    ENSG00000000460         94        102         74
    ENSG00000000938          0          0          0

#Q1.a how many genes are in this database?

``` r
nrow(counts)
```

    [1] 38694

#Q1.b How many different experiments (columns in counts or rows in
metadata) are there?

``` r
ncol(counts)
```

    [1] 8

``` r
metadata
```

              id     dex celltype     geo_id
    1 SRR1039508 control   N61311 GSM1275862
    2 SRR1039509 treated   N61311 GSM1275863
    3 SRR1039512 control  N052611 GSM1275866
    4 SRR1039513 treated  N052611 GSM1275867
    5 SRR1039516 control  N080611 GSM1275870
    6 SRR1039517 treated  N080611 GSM1275871
    7 SRR1039520 control  N061011 GSM1275874
    8 SRR1039521 treated  N061011 GSM1275875

#Q2 How many ‘control’ cell lines do we have?

``` r
sum( metadata$dex == "control")
```

    [1] 4

##4. Toy differential gene expression

To start our analysis lets calculate the mean counts for all genes in
the control experiements

1.  Extract all control columns from the counts object

2.  Calculate the mean for all rows (i.e genes) of these “control”
    columns 3-4. Do the same for “treated”

3.  compare these ‘control.mean’ and ‘treated.mean’ values.

In the lab manual we are give this chunk of code:

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
control <- metadata %>% filter(dex=="control")
control.counts <- counts %>% select(control$id) 
control.mean <- rowSums(control.counts)/4
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

#Q3 How would you make the above code in either approach more robust? is
there a function that could help here?

Instead of using `rowSums` and using a long equation, we can use
rowMeans()

``` r
#library(dplyr) don't have to load library again!
control <- metadata %>% filter(dex=="control")
control.counts <- counts %>% select(control$id) 
control.mean <- rowMeans(control.counts)
head(control.mean)
```

    ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
             900.75            0.00          520.50          339.75           97.25 
    ENSG00000000938 
               0.75 

The control columns have been extracted here:

``` r
control.inds <- metadata$dex == "control" 
control.counts <- counts [ , control.inds]
```

2.Calculate the mean for all rows (i.e genes) of these “control” columns

``` r
control.means <- rowMeans(control.counts)
```

3-4. Do the same for “treated”

#Q4 Follow the same procedure for the treated samples (i.e. calculate
the mean per gene across drug treated samples and assign to a labeled
vector called treated.mean)

``` r
treated.inds <- metadata$dex == "treated" 
```

``` r
treated.counts <- counts [ , treated.inds]
```

``` r
treated.means <- rowMeans(treated.counts)
```

#store these together for book keeping

``` r
meancounts <- data.frame(control.means, treated.means)

head(meancounts)
```

                    control.means treated.means
    ENSG00000000003        900.75        658.00
    ENSG00000000005          0.00          0.00
    ENSG00000000419        520.50        546.00
    ENSG00000000457        339.75        316.50
    ENSG00000000460         97.25         78.75
    ENSG00000000938          0.75          0.00

#Q5 (a). Create a scatter plot showing the mean of the treated samples
against the mean of the control samples. Your plot should look something
like the following.

``` r
plot(meancounts)
```

![](class12quratorr.markdown_strict_files/figure-markdown_strict/unnamed-chunk-14-1.png)

According to this plot- it seems all the points have been overplotted,
and it is highly skewed. Let ammend the data using log transformations

#Q6. Try plotting both axes on a log scale. What is the argument to
plot() that allows you to do this?

``` r
plot(meancounts, log = "xy")
```

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15032 x values <= 0 omitted
    from logarithmic plot

    Warning in xy.coords(x, y, xlabel, ylabel, log): 15281 y values <= 0 omitted
    from logarithmic plot

![](class12quratorr.markdown_strict_files/figure-markdown_strict/unnamed-chunk-15-1.png)

#Q5 (b).You could also use the ggplot2 package to make this figure
producing the plot below. What geom\_?() function would you use for this
plot?

``` r
library(ggplot2)

ggplot(meancounts, aes(x= control.means, y=treated.means))+
  geom_point(color = "pink", alpha= 0.4, size= 1)+
  labs(x= "control", y= "treated")
```

![](class12quratorr.markdown_strict_files/figure-markdown_strict/unnamed-chunk-16-1.png)

Much better! :D

We often talk metrics like “log2 fold-change”

``` r
#control/ treated

log2(10/10)
```

    [1] 0

``` r
log2(10/20)
```

    [1] -1

``` r
log2(20/10)
```

    [1] 1

``` r
log2(10/40)
```

    [1] -2

``` r
log2(40/10)
```

    [1] 2

Let’s calculate log2 fold change for our treated over control mean
counts:

``` r
meancounts$log2fc <- log2(meancounts$treated.means /
  meancounts$control.means)
```

``` r
head(meancounts)
```

                    control.means treated.means      log2fc
    ENSG00000000003        900.75        658.00 -0.45303916
    ENSG00000000005          0.00          0.00         NaN
    ENSG00000000419        520.50        546.00  0.06900279
    ENSG00000000457        339.75        316.50 -0.10226805
    ENSG00000000460         97.25         78.75 -0.30441833
    ENSG00000000938          0.75          0.00        -Inf

#Q7. What is the purpose of the arr.ind argument in the which() function
call above? Why would we then take the first column of the output and
need to call the unique() function?

``` r
zero.vals <- which(meancounts[,1:2]==0, arr.ind=TRUE)

to.rm <- unique(zero.vals[,1])
mycounts <- meancounts[-to.rm,]
head(mycounts)
```

                    control.means treated.means      log2fc
    ENSG00000000003        900.75        658.00 -0.45303916
    ENSG00000000419        520.50        546.00  0.06900279
    ENSG00000000457        339.75        316.50 -0.10226805
    ENSG00000000460         97.25         78.75 -0.30441833
    ENSG00000000971       5219.00       6687.50  0.35769358
    ENSG00000001036       2327.00       1785.75 -0.38194109

The arr.ind=TRUE argument will clause which() to return both the row and
column indices (i.e. positions) where there are TRUE values. In this
case this will tell us which genes (rows) and samples (columns) have
zero counts.

A common rule of thumb is a log2 fold change cutoff +2 and -2 to call
genes “Up regulated” or “Down regulated”.

Number of “up genes”

#Q8. Using the up.ind vector above can you determine how many up
regulated genes we have at the greater than 2 fc level?

``` r
sum(meancounts$log2fc >= +2, na.rm = T)
```

    [1] 1910

#Q9. Using the down.ind vector above can you determine how many down
regulated genes we have at the greater than 2 fc level?

Number of down genes

``` r
sum(meancounts$log2fc >= -2, na.rm = T)
```

    [1] 23046

#Q10. Do you trust these results? Why or why not?

No I do not trust these results because all our analysis has been done
based on fold change. However, fold change can be large without being
statistically significant (e.g. based on p-values). We have not done
anything yet to determine whether the differences we are seeing are
significant. These results in their current form are likely to be very
misleading.

##DESeq analysis

``` r
library(DESeq2)
```

for deseq analysis we need 3 things

-count values ‘contData\` -metadata telling us about the columns in
’countData’ (‘colData’) -design of the experiement (what do you want to
compare)

Our first function from DESeq2 will setup thr input required for
analysis by storing all these things together:

``` r
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex)
```

    converting counts to integer mode

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

``` r
dds
```

    class: DESeqDataSet 
    dim: 38694 8 
    metadata(1): version
    assays(1): counts
    rownames(38694): ENSG00000000003 ENSG00000000005 ... ENSG00000283120
      ENSG00000283123
    rowData names(0):
    colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    colData names(4): id dex celltype geo_id

The main function in DESeq2 that runs the analysis is called “DESeq()”

``` r
dds <- DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

``` r
res <- results(dds)

head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 6 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj
                    <numeric>
    ENSG00000000003  0.163035
    ENSG00000000005        NA
    ENSG00000000419  0.176032
    ENSG00000000457  0.961694
    ENSG00000000460  0.815849
    ENSG00000000938        NA

``` r
36000* 0.05
```

    [1] 1800

##Volcano summary plot

this is a common summary figure from these types of experiements and
plot the log2 fold-change vs the adjusted p-value

##Save your results

To color the points we will setup a custom color vector indicating
transcripts with large fold change and significant differences between
conditions:

``` r
# Setup our custom point color vector 
mycols <- rep("gray", nrow(res))
mycols[ abs(res$log2FoldChange) > 2 ]  <- "red" 

inds <- (res$padj < 0.01) & (abs(res$log2FoldChange) > 2 )
mycols[ inds ] <- "blue"

# Volcano plot with custom colors 
plot( res$log2FoldChange,  -log(res$padj), 
 col=mycols, ylab="-Log(P-value)", xlab="Log2(FoldChange)" )

# Cut-off lines
abline(v=c(-2,2), col="gray", lty=2)
abline(h=-log(0.1), col="gray", lty=2)
```

![](class12quratorr.markdown_strict_files/figure-markdown_strict/unnamed-chunk-32-1.png)

For even more customization you might find the EnhancedVolcano
bioconductor package useful (Note. It uses ggplot under the hood):

First we will add the more understandable gene symbol names to our full
results object res as we will use this to label the most interesting
genes in our final plot.

``` r
#BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
```

    Loading required package: ggrepel

``` r
x <- as.data.frame(res)

x$log2FoldChange <- as.numeric(x$log2FoldChange)
x$pvalue <- as.numeric(x$pvalue)
x$pvalue[x$pvalue == 0] <- 1e-300

EnhancedVolcano(x,
  lab = if ("symbol" %in% colnames(x)) x$symbol else rownames(x),
  x = "log2FoldChange",
  y = "pvalue"
)
```

    Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ℹ Please use `linewidth` instead.
    ℹ The deprecated feature was likely used in the EnhancedVolcano package.
      Please report the issue to the authors.

    Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
    ℹ Please use the `linewidth` argument instead.
    ℹ The deprecated feature was likely used in the EnhancedVolcano package.
      Please report the issue to the authors.

![](class12quratorr.markdown_strict_files/figure-markdown_strict/unnamed-chunk-33-1.png)

``` r
write.csv(res, file = "my_results.csv")
```

##Add gene annotation(second part of lab)

To help make sense of our results and communicate them to toher folks we
need to add some more annotation to our main ’ res’ object.

We will use two bioconductor packages to first map IDs to different
formats including the classic gene “symbol” gene name.

BicManager::install(“AnnotationDbi”) BicManager::install(“org.Hs.eg.db”)

``` r
library(AnnotationDbi)
```


    Attaching package: 'AnnotationDbi'

    The following object is masked from 'package:dplyr':

        select

``` r
library(org.Hs.eg.db)
```

let’s see what is ‘org.Hs.eg.db’ with the columns() function

``` r
columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

We can translate or “map” id’s between any of these 26 databases using
the ‘mapIDs()’ function

``` r
library(AnnotationDbi)
library(org.Hs.eg.db)

res$symbol <- mapIds(
  x = org.Hs.eg.db,
  keys = rownames(res),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
```

    'select()' returned 1:many mapping between keys and columns

``` r
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 7 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj      symbol
                    <numeric> <character>
    ENSG00000000003  0.163035      TSPAN6
    ENSG00000000005        NA        TNMD
    ENSG00000000419  0.176032        DPM1
    ENSG00000000457  0.961694       SCYL3
    ENSG00000000460  0.815849       FIRRM
    ENSG00000000938        NA         FGR

Add the mappings for “GENENAME” and “ENTREZID” and store as
‘res*g**e**n**e**n**a**m**e*′*a**n**d*′*r**e**s*entrez’

``` r
res$entrez <- mapIds(
  x = org.Hs.eg.db,          # database
  keys = row.names(res),     # your ENSEMBL IDs
  keytype = "ENSEMBL",       # type of your keys
  column = "ENTREZID",       # what you want back
  multiVals = "first"        # how to handle duplicates
)
```

    'select()' returned 1:many mapping between keys and columns

``` r
head(res)
```

    log2 fold change (MLE): dex treated vs control 
    Wald test p-value: dex treated vs control 
    DataFrame with 6 rows and 8 columns
                      baseMean log2FoldChange     lfcSE      stat    pvalue
                     <numeric>      <numeric> <numeric> <numeric> <numeric>
    ENSG00000000003 747.194195     -0.3507030  0.168246 -2.084470 0.0371175
    ENSG00000000005   0.000000             NA        NA        NA        NA
    ENSG00000000419 520.134160      0.2061078  0.101059  2.039475 0.0414026
    ENSG00000000457 322.664844      0.0245269  0.145145  0.168982 0.8658106
    ENSG00000000460  87.682625     -0.1471420  0.257007 -0.572521 0.5669691
    ENSG00000000938   0.319167     -1.7322890  3.493601 -0.495846 0.6200029
                         padj      symbol      entrez
                    <numeric> <character> <character>
    ENSG00000000003  0.163035      TSPAN6        7105
    ENSG00000000005        NA        TNMD       64102
    ENSG00000000419  0.176032        DPM1        8813
    ENSG00000000457  0.961694       SCYL3       57147
    ENSG00000000460  0.815849       FIRRM       55732
    ENSG00000000938        NA         FGR        2268

##Pathway analysis

There are lots of bioconductor packages to do this type of analysis. For
now let’s try one called **gage** again we need to install this if we
don’t have it already.

``` r
library(gage)
library(gageData)
library(pathview)
```

To use **gage** I need two things - a named vector offold-change valuses
for our DEGs (our genes of interest) - a set of pathways or genesets to
use for annotation.

``` r
x <- c("barry" = 5, "lisa"= 10)

x
```

    barry  lisa 
        5    10 

``` r
names(x) <- c("low","high") 

x
```

     low high 
       5   10 

``` r
foldchanges <- res$log2FoldChange
names(foldchanges) <- res$symbol

head(foldchanges)
```

         TSPAN6        TNMD        DPM1       SCYL3       FIRRM         FGR 
    -0.35070302          NA  0.20610777  0.02452695 -0.14714205 -1.73228897 

``` r
library(gage)
data(kegg.sets.hs)
foldchanges <- res$log2FoldChange
names(foldchanges) <- res$entrez
keggres <- gage(foldchanges, gsets = kegg.sets.hs)
```

In our results object:

``` r
attributes(keggres)
```

    $names
    [1] "greater" "less"    "stats"  

``` r
head(keggres$less, 5 )
```

                                                             p.geomean stat.mean
    hsa05332 Graft-versus-host disease                    0.0004250461 -3.473346
    hsa04940 Type I diabetes mellitus                     0.0017820293 -3.002352
    hsa05310 Asthma                                       0.0020045888 -3.009050
    hsa04672 Intestinal immune network for IgA production 0.0060434515 -2.560547
    hsa05330 Allograft rejection                          0.0073678825 -2.501419
                                                                 p.val      q.val
    hsa05332 Graft-versus-host disease                    0.0004250461 0.09053483
    hsa04940 Type I diabetes mellitus                     0.0017820293 0.14232581
    hsa05310 Asthma                                       0.0020045888 0.14232581
    hsa04672 Intestinal immune network for IgA production 0.0060434515 0.31387180
    hsa05330 Allograft rejection                          0.0073678825 0.31387180
                                                          set.size         exp1
    hsa05332 Graft-versus-host disease                          40 0.0004250461
    hsa04940 Type I diabetes mellitus                           42 0.0017820293
    hsa05310 Asthma                                             29 0.0020045888
    hsa04672 Intestinal immune network for IgA production       47 0.0060434515
    hsa05330 Allograft rejection                                36 0.0073678825

Let’s look at one of these pathways with our genes colored up so we can
see the overlap

``` r
pathview(pathway.id = "hsa05310", gene.data= foldchanges)
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/MacMD rendering/MacPC rendering/class 12

    Info: Writing image file hsa05310.pathview.png

Add this pathway fugure to our lab report #should be an image here

![](hsa05310.pathview.png)

Save main results

``` r
write.csv(res, file = "myresults_annotated.csv")
```
