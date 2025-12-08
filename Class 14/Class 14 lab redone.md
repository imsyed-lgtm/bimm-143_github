# Class 14: RNASeq mini project
Iman Syed A18596789

\##Background

Analysis of high-throughput biological data typically yields a list of
genes or proteins requiring further interpretation - for example the
ranked lists of differentially expressed genes we have been generating
from our RNA-seq analysis to date.

However, in many cases these ‘raw’ gene lists are challenging to
interpret due to their large size and lack of useful annotations. Hence,
our expensively assembled gene lists often fail to convey the full
degree of possible insight about the condition being studied.

Pathway analysis (also known as gene set analysis or over-representation
analysis), aims to reduce the complexity of interpreting gene lists via
mapping the listed genes to known (i.e. annotated) biological pathways,
processes and functions.

\##Data Import

Reading the `counts` and `metadata` CSV files

``` r
library(DESeq2)
```

    Loading required package: S4Vectors

    Loading required package: stats4

    Loading required package: BiocGenerics

    Loading required package: generics


    Attaching package: 'generics'

    The following objects are masked from 'package:base':

        as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
        setequal, union


    Attaching package: 'BiocGenerics'

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
        unsplit, which.max, which.min


    Attaching package: 'S4Vectors'

    The following object is masked from 'package:utils':

        findMatches

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges

    Loading required package: GenomicRanges

    Loading required package: Seqinfo

    Loading required package: SummarizedExperiment

    Loading required package: MatrixGenerics

    Loading required package: matrixStats


    Attaching package: 'MatrixGenerics'

    The following objects are masked from 'package:matrixStats':

        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars

    Loading required package: Biobase

    Welcome to Bioconductor

        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.


    Attaching package: 'Biobase'

    The following object is masked from 'package:MatrixGenerics':

        rowMedians

    The following objects are masked from 'package:matrixStats':

        anyMissing, rowMedians

``` r
metaFile <- "GSE37704_metadata.csv"
countFile <- "GSE37704_featurecounts.csv"

# Import metadata and take a peak
colData = read.csv(metaFile, row.names=1)
head(colData)
```

                  condition
    SRR493366 control_sirna
    SRR493367 control_sirna
    SRR493368 control_sirna
    SRR493369      hoxa1_kd
    SRR493370      hoxa1_kd
    SRR493371      hoxa1_kd

``` r
# Import countdata
countData = read.csv(countFile, row.names=1)
head(countData)
```

                    length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ENSG00000186092    918         0         0         0         0         0
    ENSG00000279928    718         0         0         0         0         0
    ENSG00000279457   1982        23        28        29        29        28
    ENSG00000278566    939         0         0         0         0         0
    ENSG00000273547    939         0         0         0         0         0
    ENSG00000187634   3214       124       123       205       207       212
                    SRR493371
    ENSG00000186092         0
    ENSG00000279928         0
    ENSG00000279457        46
    ENSG00000278566         0
    ENSG00000273547         0
    ENSG00000187634       258

Q1. Complete the code below to remove the troublesome first column from
countData A1:

``` r
countData <- as.matrix(countData[, -1])
head(countData)
```

                    SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
    ENSG00000186092         0         0         0         0         0         0
    ENSG00000279928         0         0         0         0         0         0
    ENSG00000279457        23        28        29        29        28        46
    ENSG00000278566         0         0         0         0         0         0
    ENSG00000273547         0         0         0         0         0         0
    ENSG00000187634       124       123       205       207       212       258

This looks better but there are lots of zero entries in there so let’s
get rid of them as we have no data for these.

Q2. Complete the code below to filter countData to exclude genes
(i.e. rows) where we have 0 read count across all samples
(i.e. columns).

Tip: What will rowSums() of countData return and how could you use it in
this context?

A2: The rowSums(countData) gives you the total count for each gene
across all samples, and genes whith all zeros will have a row sum of 0.

``` r
countData = countData[ rowSums(countData) > 0 , ]
head(countData)
```

                    SRR493366 SRR493367 SRR493368 SRR493369 SRR493370 SRR493371
    ENSG00000279457        23        28        29        29        28        46
    ENSG00000187634       124       123       205       207       212       258
    ENSG00000188976      1637      1831      2383      1226      1326      1504
    ENSG00000187961       120       153       180       236       255       357
    ENSG00000187583        24        48        65        44        48        64
    ENSG00000187642         4         9        16        14        16        16

\##DESeq analysis

Running DESeq2 Nice now lets setup the DESeqDataSet object required for
the DESeq() function and then run the DESeq pipeline. This is again
similar to our last days hands-on session.

load the packages, setup DEseq object:

``` r
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
```

    Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    design formula are characters, converting to factors

``` r
dds = DESeq(dds)
```

    estimating size factors

    estimating dispersions

    gene-wise dispersion estimates

    mean-dispersion relationship

    final dispersion estimates

    fitting model and testing

run DESeq and get results:

``` r
dds
```

    class: DESeqDataSet 
    dim: 15975 6 
    metadata(1): version
    assays(4): counts mu H cooks
    rownames(15975): ENSG00000279457 ENSG00000187634 ... ENSG00000276345
      ENSG00000271254
    rowData names(22): baseMean baseVar ... deviance maxCooks
    colnames(6): SRR493366 SRR493367 ... SRR493370 SRR493371
    colData names(2): condition sizeFactor

``` r
res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
```

Q3. Call the summary() function on your results to get a sense of how
many genes are up or down-regulated at the default 0.1 p-value cutoff.
A3: There are 4349 Up regulated genes and 4396 down reguated genes

``` r
summary(res)
```


    out of 15975 with nonzero total read count
    adjusted p-value < 0.1
    LFC > 0 (up)       : 4349, 27%
    LFC < 0 (down)     : 4396, 28%
    outliers [1]       : 0, 0%
    low counts [2]     : 1237, 7.7%
    (mean count < 0)
    [1] see 'cooksCutoff' argument of ?results
    [2] see 'independentFiltering' argument of ?results

\##Data Visualiszation Volcano plot

``` r
plot( res$log2FoldChange, -log(res$padj) )
```

![](Class-14-lab-redone_files/figure-commonmark/unnamed-chunk-9-1.png)

Q4. Improve this plot by completing the below code, which adds color and
axis labels A4:

``` r
# Make a color vector for all genes
mycols <- rep("gray", nrow(res) )

# Color red the genes with absolute fold change above 2
mycols[ abs(res$log2FoldChange) > 2 ] <- "red"

# Color blue those with adjusted p-value less than 0.01
#  and absolute fold change more than 2
inds <- (res$padj < 0.01) & (abs(res$log2FoldChange) > 2 )
mycols[ inds ] <- "blue"

plot( res$log2FoldChange, -log(res$padj), col= mycols, xlab="Log2(FoldChange)", ylab="-Log(P-value)" )
```

![](Class-14-lab-redone_files/figure-commonmark/unnamed-chunk-10-1.png)

\##Add anotation

add gene symbols

Since we mapped and counted against the Ensembl annotation, our results
only have information about Ensembl gene IDs. However, our pathway
analysis downstream will use KEGG pathways, and genes in KEGG pathways
are annotated with Entrez gene IDs. So lets add them as we did the last
day.

Q5. Use the mapIDs() function multiple times to add SYMBOL, ENTREZID and
GENENAME annotation to our results by completing the code below. A5:

``` r
library("AnnotationDbi")
library("org.Hs.eg.db")
```

``` r
columns(org.Hs.eg.db)
```

     [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
     [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
    [11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
    [16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
    [21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    [26] "UNIPROT"     

``` r
res$symbol = mapIds(org.Hs.eg.db,
                    keys = row.names(res), 
                    keytype = "ENSEMBL",
                    column = "SYMBOL",
                    multiVals = "first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$entrez = mapIds(org.Hs.eg.db,
                    keys = row.names(res),
                    keytype = "ENSEMBL",
                    column = "ENTREZID",
                    multiVals = "first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
res$name =   mapIds(org.Hs.eg.db,
                    keys = row.names(res),
                    keytype = "ENSEMBL",
                    column = "GENENAME",
                    multiVals = "first")
```

    'select()' returned 1:many mapping between keys and columns

``` r
head(res, 10)
```

    log2 fold change (MLE): condition hoxa1_kd vs control_sirna 
    Wald test p-value: condition hoxa1 kd vs control sirna 
    DataFrame with 10 rows and 9 columns
                       baseMean log2FoldChange     lfcSE       stat      pvalue
                      <numeric>      <numeric> <numeric>  <numeric>   <numeric>
    ENSG00000279457   29.913579      0.1792571 0.3248216   0.551863 5.81042e-01
    ENSG00000187634  183.229650      0.4264571 0.1402658   3.040350 2.36304e-03
    ENSG00000188976 1651.188076     -0.6927205 0.0548465 -12.630158 1.43990e-36
    ENSG00000187961  209.637938      0.7297556 0.1318599   5.534326 3.12428e-08
    ENSG00000187583   47.255123      0.0405765 0.2718928   0.149237 8.81366e-01
    ENSG00000187642   11.979750      0.5428105 0.5215598   1.040744 2.97994e-01
    ENSG00000188290  108.922128      2.0570638 0.1969053  10.446970 1.51282e-25
    ENSG00000187608  350.716868      0.2573837 0.1027266   2.505522 1.22271e-02
    ENSG00000188157 9128.439422      0.3899088 0.0467163   8.346304 7.04321e-17
    ENSG00000237330    0.158192      0.7859552 4.0804729   0.192614 8.47261e-01
                           padj      symbol      entrez                   name
                      <numeric> <character> <character>            <character>
    ENSG00000279457 6.86555e-01          NA          NA                     NA
    ENSG00000187634 5.15718e-03      SAMD11      148398 sterile alpha motif ..
    ENSG00000188976 1.76549e-35       NOC2L       26155 NOC2 like nucleolar ..
    ENSG00000187961 1.13413e-07      KLHL17      339451 kelch like family me..
    ENSG00000187583 9.19031e-01     PLEKHN1       84069 pleckstrin homology ..
    ENSG00000187642 4.03379e-01       PERM1       84808 PPARGC1 and ESRR ind..
    ENSG00000188290 1.30538e-24        HES4       57801 hes family bHLH tran..
    ENSG00000187608 2.37452e-02       ISG15        9636 ISG15 ubiquitin like..
    ENSG00000188157 4.21963e-16        AGRN      375790                  agrin
    ENSG00000237330          NA      RNF223      401934 ring finger protein ..

Q6. Finally for this section let’s reorder these results by adjusted
p-value and save them to a CSV file in your current project directory.
A6:

``` r
res = res[order(res$padj),]
write.csv(res, file="deseq_results.csv")
```

\##Section 2. Pathway Analysis Here we are going to use the gage package
for pathway analysis. Once we have a list of enriched pathways, we’re
going to use the pathview package to draw pathway diagrams, shading the
molecules in the pathway by their degree of up/down-regulation.

``` r
# Run in your R console (i.e. not your Rmarkdown doc!)
# ----BiocManager::install( c("pathview", "gage", "gageData") )

# For old vesrsions of R only (R < 3.5.0)!
#source("http://bioconductor.org/biocLite.R")
#biocLite( c("pathview", "gage", "gageData") )
```

``` r
library(pathview)
```

    ##############################################################################
    Pathview is an open source software package distributed under GNU General
    Public License version 3 (GPLv3). Details of GPLv3 is available at
    http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    formally cite the original Pathview paper (not just mention it) in publications
    or products. For details, do citation("pathview") within R.

    The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ##############################################################################

``` r
library(gage)
```

``` r
library(gageData)

data(kegg.sets.hs)
data(sigmet.idx.hs)

# Focus on signaling and metabolic pathways only
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

# Examine the first 3 pathways
head(kegg.sets.hs, 3)
```

    $`hsa00232 Caffeine metabolism`
    [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   

    $`hsa00983 Drug metabolism - other enzymes`
     [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"   "1551"  
     [9] "1553"   "1576"   "1577"   "1806"   "1807"   "1890"   "221223" "2990"  
    [17] "3251"   "3614"   "3615"   "3704"   "51733"  "54490"  "54575"  "54576" 
    [25] "54577"  "54578"  "54579"  "54600"  "54657"  "54658"  "54659"  "54963" 
    [33] "574537" "64816"  "7083"   "7084"   "7172"   "7363"   "7364"   "7365"  
    [41] "7366"   "7367"   "7371"   "7372"   "7378"   "7498"   "79799"  "83549" 
    [49] "8824"   "8833"   "9"      "978"   

    $`hsa00230 Purine metabolism`
      [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"    "10714" 
      [9] "108"    "10846"  "109"    "111"    "11128"  "11164"  "112"    "113"   
     [17] "114"    "115"    "122481" "122622" "124583" "132"    "158"    "159"   
     [25] "1633"   "171568" "1716"   "196883" "203"    "204"    "205"    "221823"
     [33] "2272"   "22978"  "23649"  "246721" "25885"  "2618"   "26289"  "270"   
     [41] "271"    "27115"  "272"    "2766"   "2977"   "2982"   "2983"   "2984"  
     [49] "2986"   "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
     [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"   "4831"  
     [65] "4832"   "4833"   "4860"   "4881"   "4882"   "4907"   "50484"  "50940" 
     [73] "51082"  "51251"  "51292"  "5136"   "5137"   "5138"   "5139"   "5140"  
     [81] "5141"   "5142"   "5143"   "5144"   "5145"   "5146"   "5147"   "5148"  
     [89] "5149"   "5150"   "5151"   "5152"   "5153"   "5158"   "5167"   "5169"  
     [97] "51728"  "5198"   "5236"   "5313"   "5315"   "53343"  "54107"  "5422"  
    [105] "5424"   "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
    [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"   "5441"  
    [121] "5471"   "548644" "55276"  "5557"   "5558"   "55703"  "55811"  "55821" 
    [129] "5631"   "5634"   "56655"  "56953"  "56985"  "57804"  "58497"  "6240"  
    [137] "6241"   "64425"  "646625" "654364" "661"    "7498"   "8382"   "84172" 
    [145] "84265"  "84284"  "84618"  "8622"   "8654"   "87178"  "8833"   "9060"  
    [153] "9061"   "93034"  "953"    "9533"   "954"    "955"    "956"    "957"   
    [161] "9583"   "9615"  

The main gage() function requires a named vector of fold changes, where
the names of the values are the Entrez gene IDs.

Note that we used the mapIDs() function above to obtain Entrez gene IDs
(stored in
res$entrez) and we have the fold change results from DESeq2 analysis (stored in res$log2FoldChange).

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

         1266     54855      1465      2034      2150      6659 
    -2.422719  3.201955 -2.313738 -1.888019  3.344508  2.392288 

Now, let’s run the gage pathway analysis.

``` r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

Now lets look at the object returned from gage().

``` r
attributes(keggres)
```

    $names
    [1] "greater" "less"    "stats"  

``` r
head(keggres$less)
```

                                             p.geomean stat.mean        p.val
    hsa04110 Cell cycle                   8.995727e-06 -4.378644 8.995727e-06
    hsa03030 DNA replication              9.424076e-05 -3.951803 9.424076e-05
    hsa03013 RNA transport                1.375901e-03 -3.028500 1.375901e-03
    hsa03440 Homologous recombination     3.066756e-03 -2.852899 3.066756e-03
    hsa04114 Oocyte meiosis               3.784520e-03 -2.698128 3.784520e-03
    hsa00010 Glycolysis / Gluconeogenesis 8.961413e-03 -2.405398 8.961413e-03
                                                q.val set.size         exp1
    hsa04110 Cell cycle                   0.001448312      121 8.995727e-06
    hsa03030 DNA replication              0.007586381       36 9.424076e-05
    hsa03013 RNA transport                0.073840037      144 1.375901e-03
    hsa03440 Homologous recombination     0.121861535       28 3.066756e-03
    hsa04114 Oocyte meiosis               0.121861535      102 3.784520e-03
    hsa00010 Glycolysis / Gluconeogenesis 0.212222694       53 8.961413e-03

Now, let’s try out the pathview() function from the pathview package to
make a pathway plot with our RNA-Seq expression results shown in color.
To begin with lets manually supply a pathway.id (namely the first part
of the “hsa04110 Cell cycle”) that we could see from the print out
above.

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04110.pathview.png

![](hsa04110.pathview.png) This downloads the pathway figure data from
KEGG and adds our results to it. Here is the default low resolution
raster PNG output from the pathview() call above:

``` r
# A different PDF based output of the same data
pathview(gene.data=foldchanges, pathway.id="hsa04110", kegg.native=FALSE)
```

    'select()' returned 1:1 mapping between keys and columns

    Warning: reconcile groups sharing member nodes!

         [,1] [,2] 
    [1,] "9"  "300"
    [2,] "9"  "306"

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04110.pathview.pdf

Now, let’s process our results a bit more to automagicaly pull out the
top 5 upregulated pathways, then further process that just to get the
pathway IDs needed by the pathview() function. We’ll use these KEGG
pathway IDs for pathview plotting below.

``` r
## Focus on top 5 upregulated pathways here for demo purposes only
keggrespathways <- rownames(keggres$greater)[1:5]

# Extract the 8 character long IDs part of each string
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
```

    [1] "hsa04640" "hsa04630" "hsa00140" "hsa04142" "hsa04330"

Finally, lets pass these IDs in keggresids to the pathview() function to
draw plots for all the top 5 pathways.

``` r
pathview(gene.data=foldchanges, pathway.id=keggresids, species="hsa")
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04640.pathview.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04630.pathview.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa00140.pathview.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04142.pathview.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04330.pathview.png

![](hsa04640.pathview.png) ![](hsa00140.pathview.png)
![](hsa04142.pathview.png) ![](hsa04330.pathview.png)
![](hsa04630.pathview.png)

Q7. Can you do the same procedure as above to plot the pathview figures
for the top 5 down-reguled pathways? A7: (use the same code for up
regulated pathways, but just change to down regulated)

``` r
keggres_down <- rownames(keggres$less)[1:5]
keggresids_down <- substr(keggres_down, 1, 8)  # first 8 characters = KEGG IDs

# 3️⃣ Generate PNGs in the current working directory
lapply(keggresids_down, function(pid) {
  pathview(
    gene.data = foldchanges,
    pathway.id = pid,
    species = "hsa",
    kegg.native = TRUE,
    out.suffix = pid       # ensures unique filenames
    # out.dir is omitted, so PNGs go to current working directory
  )
})
```

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04110.hsa04110.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa03030.hsa03030.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa03013.hsa03013.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa03440.hsa03440.png

    'select()' returned 1:1 mapping between keys and columns

    Info: Working in directory /Users/imsyed/Documents/bimm 143 class 14

    Info: Writing image file hsa04114.hsa04114.png

    [[1]]
    [[1]]$plot.data.gene
        kegg.names   labels
    4         1029   CDKN2A
    5        51343     FZR1
    6         4171     MCM2
    7         4998     ORC1
    8          996    CDC27
    9          996    CDC27
    10        6500     SKP1
    11        6500     SKP1
    24         983     CDK1
    25         701    BUB1B
    26        5111     PCNA
    27        5347     PLK1
    28         472      ATM
    29       10926     DBF4
    30        4085   MAD2L1
    31         699     BUB1
    32        9184     BUB3
    33        8555   CDC14B
    34        7529    YWHAB
    35        2810      SFN
    36        1111    CHEK1
    37        1026   CDKN1A
    38        5591    PRKDC
    39        4193     MDM2
    40        1387   CREBBP
    41        6502     SKP2
    42        9088   PKMYT1
    43        7465     WEE1
    44        9232    PTTG1
    45        9700    ESPL1
    46        8243    SMC1A
    48        6502     SKP2
    49        5925      RB1
    50        1647  GADD45A
    51        5925      RB1
    52        7157     TP53
    53        1027   CDKN1B
    54        1030   CDKN2B
    55        7040    TGFB1
    56        4089    SMAD4
    57        4087    SMAD2
    58        8317     CDC7
    60         991    CDC20
    61         994   CDC25B
    62        8318    CDC45
    63         990     CDC6
    64         993   CDC25A
    65        2932    GSK3B
    66         983     CDK1
    67         891    CCNB1
    68         890    CCNA2
    69        1022     CDK7
    70        1017     CDK2
    71        1017     CDK2
    72        1019     CDK4
    73         902     CCNH
    74         890    CCNA2
    75         898    CCNE1
    76         595    CCND1
    79        8379   MAD1L1
    82        7272      TTK
    175       1031   CDKN2C
    176       1032   CDKN2D
    177       1029   CDKN2A
    178       9126     SMC3
    179      10274    STAG1
    180       5885    RAD21
    186       5933     RBL1
    187       1874     E2F4
    188       4609      MYC
    189       7709   ZBTB17
    196         25     ABL1
    197       7027    TFDP1
    198       1869     E2F1
    199       3065    HDAC1
    200       5925      RB1
    201       5933     RBL1
    214       7027    TFDP1
    215       1874     E2F4
    218       7027    TFDP1
    225      10403    NDC80
    226       9212    AURKB
    227      57082     KNL1
    228       9587 MAD2L1BP
    230       9319   TRIP13
    237       8243    SMC1A
    238       9126     SMC3
    239        546     ATRX
    240      25836    NIPBL
    241      23383     MAU2
    242      81620     CDT1
    243      27085     MTBP
    244      90381    TICRR
    246     114799    ESCO1
    248       5885    RAD21
    249      10274    STAG1
    251       1663    DDX11
    252       8243    SMC1A
    253       9126     SMC3
    254      23047    PDS5B
    255      23063     WAPL
    256       5885    RAD21
    257      10274    STAG1
    260     113130    CDCA5
    261      55869    HDAC8
    264        983     CDK1
    266       5347     PLK1
    276     151648     SGO1
    277       5515   PPP2CA
    283       1026   CDKN1A
    284      26271    FBXO5
    289       7157     TP53
                                                                               all.mapped
    4                                                                                1029
    5                                                                               51343
    6                                                       4171,4172,4173,4174,4175,4176
    7                                                     4998,4999,5000,5001,23594,23595
    8   996,8697,8881,10393,25847,25906,29882,29945,51433,51434,51529,64682,119504,246184
    9   996,8697,8881,10393,25847,25906,29882,29945,51433,51434,51529,64682,119504,246184
    10                                                                     6500,8454,9978
    11                                                                     6500,8454,9978
    24                                                                                983
    25                                                                                701
    26                                                                               5111
    27                                                                               5347
    28                                                                            472,545
    29                                                                        10926,80174
    30                                                                         4085,10459
    31                                                                                699
    32                                                                               9184
    33                                                                          8555,8556
    34                                                     7529,7531,7532,7533,7534,10971
    35                                                                                   
    36                                                                         1111,11200
    37                                                                               1026
    38                                                                               5591
    39                                                                               4193
    40                                                                          1387,2033
    41                                                                               6502
    42                                                                               9088
    43                                                                        7465,494551
    44                                                                               9232
    45                                                                               9700
    46                                                                               8243
    48                                                                               6502
    49                                                                               5925
    50                                                                    1647,4616,10912
    51                                                                               5925
    52                                                                               7157
    53                                                                          1027,1028
    54                                                                               1030
    55                                                                     7040,7042,7043
    56                                                                               4089
    57                                                                          4087,4088
    58                                                                               8317
    60                                                                                991
    61                                                                            994,995
    62                                                                               8318
    63                                                                                990
    64                                                                                993
    65                                                                               2932
    66                                                                                983
    67                                                                     891,9133,85417
    68                                                                           890,8900
    69                                                                               1022
    70                                                                               1017
    71                                                                               1017
    72                                                                          1019,1021
    73                                                                                902
    74                                                                           890,8900
    75                                                                           898,9134
    76                                                                        595,894,896
    79                                                                               8379
    82                                                                               7272
    175                                                                              1031
    176                                                                              1032
    177                                                                              1029
    178                                                                              9126
    179                                                                       10274,10735
    180                                                                              5885
    186                                                                              5933
    187                                                                         1874,1875
    188                                                                              4609
    189                                                                              7709
    196                                                                                25
    197                                                                         7027,7029
    198                                                                    1869,1870,1871
    199                                                                         3065,3066
    200                                                                              5925
    201                                                                         5933,5934
    214                                                                         7027,7029
    215                                                                         1874,1875
    218                                                                         7027,7029
    225                                                                             10403
    226                                                                              9212
    227                                                                             57082
    228                                                                              9587
    230                                                                              9319
    237                                                                              8243
    238                                                                              9126
    239                                                                               546
    240                                                                             25836
    241                                                                             23383
    242                                                                             81620
    243                                                                             27085
    244                                                                             90381
    246                                                                     114799,157570
    248                                                                              5885
    249                                                                       10274,10735
    251                                                                              1663
    252                                                                              8243
    253                                                                              9126
    254                                                                       23047,23244
    255                                                                             23063
    256                                                                              5885
    257                                                                       10274,10735
    260                                                                            113130
    261                                                                             55869
    264                                                                               983
    266                                                                              5347
    276                                                                            151648
    277                                      5515,5516,5518,5519,5525,5526,5527,5528,5529
    283                                                                              1026
    284                                                                             26271
    289                                                                              7157
        type    x   y width height     mol.data mol.col
    4   gene  532 218    46     17  0.119338916 #BEBEBE
    5   gene  981 630    46     17 -0.367315682 #8FCE8F
    6   gene  553 681    46     17 -9.133771506 #00FF00
    7   gene  494 681    46     17 -5.080008665 #00FF00
    8   gene  981 392    46     17 -3.115191980 #00FF00
    9   gene  981 613    46     17 -3.115191980 #00FF00
    10  gene  188 613    46     17 -0.732831812 #30EF30
    11  gene  432 285    46     17 -0.732831812 #30EF30
    24  gene  780 562    46     17 -2.726227119 #00FF00
    25  gene  873 392    46     17 -3.085447019 #00FF00
    26  gene  559 453    46     17  0.260140443 #CE8F8F
    27  gene  862 562    46     17 -3.329237417 #00FF00
    28  gene  696 257    46     17  0.564593211 #DF5F5F
    29  gene  532 748    46     17 -2.963221676 #00FF00
    30  gene  873 375    46     17 -1.954222022 #00FF00
    31  gene  808 409    46     17 -3.032368213 #00FF00
    32  gene  873 409    46     17 -0.361252978 #8FCE8F
    33  gene 1009 672    46     17 -0.294882520 #8FCE8F
    34  gene  850 446    46     17 -2.665753282 #00FF00
    35  gene  638 407    46     17           NA #FFFFFF
    36  gene  696 393    46     17 -2.053242791 #00FF00
    37  gene  459 407    46     17  1.832257777 #FF0000
    38  gene  640 257    46     17 -0.574579028 #5FDF5F
    39  gene  532 273    46     17 -0.417188973 #5FDF5F
    40  gene  590 244    46     17  0.208942443 #CE8F8F
    41  gene  432 302    46     17 -1.738103676 #00FF00
    42  gene  763 622    46     17 -2.703573241 #00FF00
    43  gene  704 622    46     17 -0.382219288 #8FCE8F
    44  gene  981 291    46     17 -2.766072229 #00FF00
    45  gene  981 204    46     17 -2.745492001 #00FF00
    46  gene  957  79    46     17 -0.646639657 #30EF30
    48  gene  188 630    46     17 -1.738103676 #00FF00
    49  gene  637 622    46     17  0.356540702 #CE8F8F
    50  gene  559 407    46     17 -1.198361637 #00FF00
    51  gene  504 337    46     17  0.356540702 #CE8F8F
    52  gene  590 337    46     17 -0.335083507 #8FCE8F
    53  gene  407 407    46     17  1.307620792 #FF0000
    54  gene  256 407    46     17  2.179829794 #FF0000
    55  gene  327 198    46     17  2.126655915 #FF0000
    56  gene  327 258    46     17  0.092515835 #BEBEBE
    57  gene  327 241    46     17  0.207771865 #CE8F8F
    58  gene  532 731    46     17 -1.843537902 #00FF00
    60  gene  981 409    46     17 -3.233403735 #00FF00
    61  gene  830 496    46     17 -3.158955125 #00FF00
    62  gene  553 638    46     17 -2.474892600 #00FF00
    63  gene  494 647    46     17 -2.139277489 #00FF00
    64  gene  614 496    46     17 -1.539576694 #00FF00
    65  gene  170 261    46     17  0.533894665 #DF5F5F
    66  gene  689 562    46     17 -2.726227119 #00FF00
    67  gene  780 545    46     17 -6.006663761 #00FF00
    68  gene  689 545    46     17 -0.845750930 #00FF00
    69  gene  552 562    46     17  0.061254714 #BEBEBE
    70  gene  457 562    46     17  0.834646948 #FF0000
    71  gene  335 562    46     17  0.834646948 #FF0000
    72  gene  227 562    46     17  0.713179703 #EF3030
    73  gene  552 545    46     17  0.068366475 #BEBEBE
    74  gene  457 545    46     17 -0.845750930 #00FF00
    75  gene  335 545    46     17 -1.619361532 #00FF00
    76  gene  227 545    46     17 -0.410015021 #5FDF5F
    79  gene  874 325    46     17  0.007601264 #BEBEBE
    82  gene  874 257    46     17 -3.159407897 #00FF00
    175 gene  304 407    46     17 -2.757140189 #00FF00
    176 gene  352 407    46     17  0.149572876 #BEBEBE
    177 gene  208 407    46     17  0.119338916 #BEBEBE
    178 gene 1003  79    46     17 -0.907010638 #00FF00
    179 gene  980 113    46     17  0.797546572 #EF3030
    180 gene  980  96    46     17 -0.420884017 #5FDF5F
    186 gene  243 224    46     17 -1.808482919 #00FF00
    187 gene  243 241    46     17 -0.183213204 #BEBEBE
    188 gene  243 313    46     17  0.157441084 #BEBEBE
    189 gene  243 345    46     17  0.097634890 #BEBEBE
    196 gene  413 620    46     17 -0.244215556 #8FCE8F
    197 gene  353 688    46     17  0.357053597 #CE8F8F
    198 gene  353 671    46     17 -3.253523136 #00FF00
    199 gene  413 644    46     17  0.500255529 #DF5F5F
    200 gene  353 630    46     17  0.356540702 #CE8F8F
    201 gene  298 630    46     17 -1.490970660 #00FF00
    214 gene  298 688    46     17  0.357053597 #CE8F8F
    215 gene  298 671    46     17 -0.183213204 #BEBEBE
    218 gene  243 258    46     17  0.357053597 #CE8F8F
    225 gene  808 257    46     17 -2.864997646 #00FF00
    226 gene  808 129    46     17 -3.064158113 #00FF00
    227 gene  808 364    46     17 -3.002387278 #00FF00
    228 gene  930 342    46     17 -0.027540435 #BEBEBE
    230 gene  930 325    46     17 -2.537246465 #00FF00
    237 gene  219  79    46     17 -0.646639657 #30EF30
    238 gene  265  79    46     17 -0.907010638 #00FF00
    239 gene  269  54    46     17  0.403085300 #DF5F5F
    240 gene  327  79    46     17  0.510547216 #DF5F5F
    241 gene  327  96    46     17  0.288305186 #CE8F8F
    242 gene  523 664    46     17 -1.517394041 #00FF00
    243 gene  553 621    46     17 -1.173272386 #00FF00
    244 gene  553 604    46     17 -1.909648862 #00FF00
    246 gene  613  54    46     17 -3.041352707 #00FF00
    248 gene  243  96    46     17 -0.420884017 #5FDF5F
    249 gene  243 113    46     17  0.797546572 #EF3030
    251 gene  320  54    46     17 -2.077638964 #00FF00
    252 gene  494  79    46     17 -0.646639657 #30EF30
    253 gene  540  79    46     17 -0.907010638 #00FF00
    254 gene  590 104    46     17  0.298149305 #CE8F8F
    255 gene  590 121    46     17 -0.058150017 #BEBEBE
    256 gene  518  96    46     17 -0.420884017 #5FDF5F
    257 gene  518 113    46     17  0.797546572 #EF3030
    260 gene  662 104    46     17 -2.272633833 #00FF00
    261 gene  540  54    46     17 -0.153951685 #BEBEBE
    264 gene  808 104    46     17 -2.726227119 #00FF00
    266 gene  808  79    46     17 -3.329237417 #00FF00
    276 gene  913 104    46     17 -2.654160384 #00FF00
    277 gene  913 121    46     17 -2.779166545 #00FF00
    283 gene  953 714    46     17  1.832257777 #FF0000
    284 gene  953 672    46     17 -1.413543648 #00FF00
    289 gene  954 756    46     17 -0.335083507 #8FCE8F

    [[1]]$plot.data.cpd
    NULL


    [[2]]
    [[2]]$plot.data.gene
        kegg.names   labels all.mapped type    x    y width height    mol.data
    4         3978     LIG1       3978 gene 1158 1231    46     17 -0.84011948
    5         2237     FEN1       2237 gene 1081 1231    46     17 -1.20475314
    6         1763     DNA2       1763 gene 1003 1230    46     17 -0.30606476
    7        10535 RNASEH2A      10535 gene 1067 1174    46     17 -1.71320965
    8        79621 RNASEH2B      79621 gene 1113 1174    46     17 -0.50604678
    9         5982     RFC2  5982,5984 gene 1113 1116    46     17 -1.77251145
    10       84153 RNASEH2C      84153 gene 1159 1174    46     17 -0.43587616
    11        5983     RFC3  5983,5985 gene 1159 1116    46     17 -1.99908008
    12        5111     PCNA       5111 gene 1003 1116    46     17  0.26014044
    13        6119     RPA3       6119 gene 1159 1065    46     17 -1.00519942
    14        4174     MCM5       4174 gene 1068 1048    46     17 -2.46315174
    15       57804    POLD4      57804 gene 1159  918    46     17  0.88890207
    16       10714    POLD3      10714 gene 1113  918    46     17 -0.52707926
    17       56655    POLE4      56655 gene 1159  976    46     17  0.03962253
    18       54107    POLE3      54107 gene 1113  976    46     17 -0.37024407
    19        6742    SSBP1       6742 gene 1158  309    46     17 -0.18171090
    31      246243  RNASEH1     246243 gene 1003  366    46     17 -0.02864004
    33        4171     MCM2       4171 gene 1022 1031    46     17 -1.80843478
    34        4172     MCM3       4172 gene 1068 1031    46     17 -1.11789555
    35        4173     MCM4       4173 gene 1022 1048    46     17 -0.67137753
    36        5422    POLA1       5422 gene 1021  862    46     17 -1.27047988
    37       23649    POLA2      23649 gene 1067  862    46     17 -1.17429437
    38        5557    PRIM1       5557 gene 1113  862    46     17 -2.69232304
    39        5558    PRIM2       5558 gene 1159  862    46     17 -1.14726049
    40        5424    POLD1       5424 gene 1021  918    46     17 -1.10108748
    41        5425    POLD2       5425 gene 1067  918    46     17 -0.68208724
    42        5426     POLE       5426 gene 1021  976    46     17 -0.95033100
    43        5427    POLE2       5427 gene 1067  976    46     17 -2.64019397
    44        4175     MCM6       4175 gene 1022 1065    46     17 -2.10497621
    45        4176     MCM7       4176 gene 1068 1065    46     17 -0.96793569
    46        6117     RPA1       6117 gene 1159 1031    46     17 -0.77769183
    47        6118     RPA2 6118,29935 gene 1159 1048    46     17  0.91159289
    48        5981     RFC1       5981 gene 1067 1116    46     17 -0.86785058
    96        6117     RPA1       6117 gene 1158  609    46     17 -0.77769183
    97        5111     PCNA       5111 gene 1003  679    46     17  0.26014044
    99      246243  RNASEH1     246243 gene 1158  679    46     17 -0.02864004
    100       3978     LIG1       3978 gene 1158  753    46     17 -0.84011948
    101       2237     FEN1       2237 gene 1081  753    46     18 -1.20475314
    103       1763     DNA2       1763 gene 1003  753    46     17 -0.30606476
    105     246243  RNASEH1     246243 gene 1003 1174    46     17 -0.02864004
        mol.col
    4   #00FF00
    5   #00FF00
    6   #8FCE8F
    7   #00FF00
    8   #5FDF5F
    9   #00FF00
    10  #5FDF5F
    11  #00FF00
    12  #CE8F8F
    13  #00FF00
    14  #00FF00
    15  #FF0000
    16  #5FDF5F
    17  #BEBEBE
    18  #8FCE8F
    19  #BEBEBE
    31  #BEBEBE
    33  #00FF00
    34  #00FF00
    35  #30EF30
    36  #00FF00
    37  #00FF00
    38  #00FF00
    39  #00FF00
    40  #00FF00
    41  #30EF30
    42  #00FF00
    43  #00FF00
    44  #00FF00
    45  #00FF00
    46  #30EF30
    47  #FF0000
    48  #00FF00
    96  #30EF30
    97  #CE8F8F
    99  #BEBEBE
    100 #00FF00
    101 #00FF00
    103 #8FCE8F
    105 #BEBEBE

    [[2]]$plot.data.cpd
    NULL


    [[3]]
    [[3]]$plot.data.gene
        kegg.names  labels
    442       1434   CSE1L
    443       5901     RAN
    446       4686   NCBP1
    447      51808    PHAX
    530       9818   NUP58
    531      23165  NUP205
    532      23511  NUP188
    533      23636   NUP62
    534       7175     TPR
    535       9972  NUP153
    536       9688   NUP93
    537      53371   NUP54
    538      10762   NUP50
    539       9631  NUP155
    540      59343   SENP2
    541       4928   NUP98
    542      79023   NUP37
    543       6396   SEC13
    544      23279  NUP160
    545      57122  NUP107
    546      79902   NUP85
    547      81929   SEH1L
    548       4928   NUP98
    549       8480    RAE1
    550     348995   NUP43
    551      55746  NUP133
    552      23225  NUP210
    553       9883  POM121
    554      55706    NDC1
    555       8021  NUP214
    556       4927   NUP88
    557       5903  RANBP2
    558       7329   UBE2I
    559       5905 RANGAP1
    560       6612   SUMO3
    561      11097   NUP42
    562       8086    AAAS
    570       1915  EEF1A1
    626     129401   NUP35
    650      11269  DDX19B
    651       2733    GLE1
    677      25909  AHCTF1
    679      55161  TMEM33
    685      51068    NMD3
    689       9939   RBM8A
    690       4116   MAGOH
    691       9775  EIF4A3
    692      22794   CASC3
    693      22985   ACIN1
    694      10921   RNPS1
    695      10284   SAP18
    696      10250   SRRM1
    697      10482    NXF1
    698       5411     PNN
    699      10189  ALYREF
    700      29107    NXT1
    701       5976    UPF1
    702       7919  DDX39B
    703       9984   THOC1
    704      57187   THOC2
    705      84321   THOC3
    706       8563   THOC5
    707      79228   THOC6
    708      80145   THOC7
    711      26019    UPF2
    712      65109   UPF3B
    715      84305    PYM1
    720       3837   KPNB1
    721      10073   SNUPN
    731       3836   KPNA1
                                                           all.mapped type   x   y
    442                 1434,7514,11260,23039,23214,57510,64328,64901 gene 495 798
    443                                                          5901 gene 541 798
    446                                                    4686,22916 gene 641 822
    447                                                         51808 gene 595 822
    530                                                          9818 gene 111 719
    531                                                         23165 gene 211 719
    532                                                         23511 gene 257 719
    533                                                         23636 gene  65 719
    534                                                          7175 gene  65 818
    535                                                          9972 gene 176 819
    536                                                          9688 gene 349 719
    537                                                         53371 gene 157 719
    538                                                         10762 gene 121 818
    539                                                          9631 gene 303 719
    540                                                         59343 gene 222 819
    541                                                          4928 gene 311 658
    542                                                         79023 gene 458 658
    543                                                          6396 gene 156 658
    544                                                         23279 gene  65 658
    545                                                         57122 gene 202 658
    546                                                         79902 gene 111 658
    547                                                         81929 gene 366 658
    548                                                          4928 gene 310 604
    549                                                          8480 gene 264 604
    550                                                        348995 gene 412 658
    551                                                         55746 gene 248 658
    552                                                   23225,91181 gene 116 771
    553                                          9883,94026,100101267 gene 166 771
    554                                                         55706 gene  65 771
    555                                                          8021 gene 361 604
    556                                                          4927 gene 407 604
    557                                                          5903 gene 468 603
    558                                                          7329 gene 560 603
    559                                                          5905 gene 514 603
    560                                         6612,6613,7341,387082 gene 606 603
    561                                                         11097 gene 114 604
    562                                                          8086 gene  64 604
    570                                                     1915,1917 gene 595 798
    626                                                        129401 gene 394 719
    650                                                   11269,55308 gene 214 604
    651                                                          2733 gene 164 604
    677                                                         25909 gene 514 658
    679                                                         55161 gene 328 771
    685                                                         51068 gene 595 845
    689                                                          9939 gene 736 640
    690                                                    4116,55110 gene 782 640
    691                                                          9775 gene 874 640
    692                                                         22794 gene 828 640
    693                                                         22985 gene 737 687
    694                                                         10921 gene 828 687
    695                                                         10284 gene 783 687
    696                                                         10250 gene 887 762
    697                                                         10482 gene 737 762
    698                                                          5411 gene 874 687
    699                                                         10189 gene 920 687
    700                                                   29107,55916 gene 783 762
    701                                                          5976 gene 737 736
    702                                                          7919 gene 835 762
    703                                                          9984 gene 737 843
    704                                                         57187 gene 783 843
    705                                                         84321 gene 969 843
    706                                                          8563 gene 828 843
    707                                                         79228 gene 874 843
    708                                                         80145 gene 920 843
    711                                                         26019 gene 786 736
    712                                                   65109,65110 gene 835 736
    715                                                         84305 gene 939 762
    720 3837,3842,3843,9670,10526,10527,23534,30000,51194,55705,79711 gene 541 758
    721                                                         10073 gene 595 758
    731                         3836,3838,3839,3840,3841,23633,402569 gene 495 758
        width height     mol.data mol.col
    442    46     17 -5.507105307 #00FF00
    443    46     17 -0.542156931 #5FDF5F
    446    46     17  0.389826940 #CE8F8F
    447    46     17 -0.228583301 #8FCE8F
    530    46     17  0.102527952 #BEBEBE
    531    46     17 -0.765740818 #30EF30
    532    46     17 -0.264777916 #8FCE8F
    533    46     17 -0.751648232 #30EF30
    534    46     17  0.413781321 #DF5F5F
    535    46     17 -0.016390011 #BEBEBE
    536    46     17 -0.407666536 #5FDF5F
    537    46     17 -0.294131066 #8FCE8F
    538    46     17  0.461161581 #DF5F5F
    539    46     17 -0.542487363 #5FDF5F
    540    46     17 -0.439385544 #5FDF5F
    541    46     17 -0.051371108 #BEBEBE
    542    46     17 -0.592509552 #5FDF5F
    543    46     17 -0.433377920 #5FDF5F
    544    46     17  0.057617485 #BEBEBE
    545    46     17 -0.341975989 #8FCE8F
    546    46     17 -0.355402186 #8FCE8F
    547    46     17 -1.097482803 #00FF00
    548    46     17 -0.051371108 #BEBEBE
    549    46     17  0.032152967 #BEBEBE
    550    46     17 -0.708047032 #30EF30
    551    46     17  0.112137249 #BEBEBE
    552    46     17 -1.169844000 #00FF00
    553    46     17 -2.064486416 #00FF00
    554    46     17 -1.664759901 #00FF00
    555    46     17 -0.322067160 #8FCE8F
    556    46     17 -0.411414422 #5FDF5F
    557    46     17  0.260913375 #CE8F8F
    558    46     17  0.262428839 #CE8F8F
    559    46     17 -0.522909044 #5FDF5F
    560    46     17 -3.004894688 #00FF00
    561    46     17 -0.016752472 #BEBEBE
    562    46     17 -0.297953817 #8FCE8F
    570    46     17 -0.076555439 #BEBEBE
    626    46     17 -0.562388659 #5FDF5F
    650    46     17 -0.234704576 #8FCE8F
    651    46     17 -0.292745742 #8FCE8F
    677    46     17 -0.081250008 #BEBEBE
    679    46     17  0.612185021 #EF3030
    685    46     17  0.233420696 #CE8F8F
    689    46     17  0.179459515 #BEBEBE
    690    46     17 -0.844205700 #00FF00
    691    46     17 -0.078236758 #BEBEBE
    692    46     17 -0.263613986 #8FCE8F
    693    46     17  0.038279575 #BEBEBE
    694    46     17  0.236091030 #CE8F8F
    695    46     17 -0.034400880 #BEBEBE
    696    46     17 -0.517840234 #5FDF5F
    697    46     17  0.203820686 #CE8F8F
    698    46     17 -0.218475108 #8FCE8F
    699    46     17 -0.316664332 #8FCE8F
    700    46     17  0.194914502 #BEBEBE
    701    46     17  0.051029751 #BEBEBE
    702    46     17  0.488894616 #DF5F5F
    703    46     17 -0.877754507 #00FF00
    704    46     17  0.150760859 #BEBEBE
    705    46     17 -1.200370576 #00FF00
    706    46     17  0.118070499 #BEBEBE
    707    46     17  0.085989957 #BEBEBE
    708    46     17 -0.201186617 #8FCE8F
    711    46     17 -0.006487418 #BEBEBE
    712    46     17  0.613192722 #EF3030
    715    46     17  0.278869504 #CE8F8F
    720    46     17 -4.288326490 #00FF00
    721    46     17 -0.704982949 #30EF30
    731    46     17 -4.514665321 #00FF00

    [[3]]$plot.data.cpd
    NULL


    [[4]]
    [[4]]$plot.data.gene
        kegg.names   labels            all.mapped type    x   y width height
    65        7156    TOP3A             7156,8940 gene  636 756    46     17
    66         641      BLM                   641 gene  905 885    46     17
    67        7979     SEM1                  7979 gene 1087 351    46     17
    68        5893    RAD52                  5893 gene 1015 373    46     17
    69        4361    MRE11                  4361 gene  998 219    46     17
    70       10111    RAD50                 10111 gene  998 202    46     17
    71        7517    XRCC3                  7517 gene  646 358    46     17
    72        5889   RAD51C                  5889 gene  646 375    46     17
    73        7516    XRCC2                  7516 gene  747 392    46     17
    74        7516    XRCC2                  7516 gene  593 375    46     17
    75        5892   RAD51D                  5892 gene  747 375    46     17
    76        5892   RAD51D                  5892 gene  593 358    46     17
    77        5889   RAD51C                  5889 gene  701 375    46     17
    78        5890   RAD51B                  5890 gene  701 358    46     17
    79        5889   RAD51C                  5889 gene  593 326    46     17
    95       80198    MUS81                 80198 gene  688 739    46     17
    99        4361    MRE11                  4361 gene  902 219    46     17
    100      10111    RAD50                 10111 gene  902 202    46     17
    101       5890   RAD51B                  5890 gene  593 309    46     17
    102       4683      NBN                  4683 gene  998 236    46     17
    103       6117     RPA1  6117,6118,6119,29935 gene 1015 329    46     17
    104       5888    RAD51                  5888 gene 1015 351    46     17
    110     146956     EME1                146956 gene  688 756    46     17
    112        675    BRCA2                   675 gene 1133 351    46     17
    114       8438   RAD54L            8438,25788 gene 1012 493    46     17
    115       5424    POLD1 5424,5425,10714,57804 gene  917 668    46     17
    116        641      BLM                   641 gene  636 739    46     17
    120       6742    SSBP1                  6742 gene  253 202    46     17
    321      50511    SYCP3                 50511 gene 1214 351    46     17
    324        672    BRCA1                   672 gene 1133 236    46     17
    325       5932    RBBP8                  5932 gene 1087 219    46     17
    326      79728    PALB2                 79728 gene 1133 334    46     17
    327      51720    UIMC1                 51720 gene 1159 293    46     17
    328      84142 ABRAXAS1                 84142 gene 1159 276    46     17
    329      83990    BRIP1                 83990 gene 1179 236    46     17
    330        580    BARD1                   580 gene 1087 236    46     17
    331      11073   TOPBP1                 11073 gene 1179 208    46     17
    332        472      ATM                   472 gene 1133 117    46     17
    335       9577   BABAM2                  9577 gene 1214 276    46     17
    336      29086   BABAM1                 29086 gene 1214 259    46     17
    337      79184    BRCC3                 79184 gene 1214 293    46     17
           mol.data mol.col
    65  -0.61651729 #30EF30
    66  -1.64085319 #00FF00
    67   0.18551461 #BEBEBE
    68  -0.44503968 #5FDF5F
    69  -0.82731603 #00FF00
    70   0.33035621 #CE8F8F
    71  -0.47755042 #5FDF5F
    72  -0.46099838 #5FDF5F
    73  -2.10623179 #00FF00
    74  -2.10623179 #00FF00
    75  -0.64045533 #30EF30
    76  -0.64045533 #30EF30
    77  -0.46099838 #5FDF5F
    78  -0.66137734 #30EF30
    79  -0.46099838 #5FDF5F
    95   0.04536843 #BEBEBE
    99  -0.82731603 #00FF00
    100  0.33035621 #CE8F8F
    101 -0.66137734 #30EF30
    102  0.02265732 #BEBEBE
    103 -0.87129836 #00FF00
    104 -2.25179994 #00FF00
    110 -2.39887847 #00FF00
    112 -2.12177537 #00FF00
    114 -4.66279482 #00FF00
    115 -1.42135191 #00FF00
    116 -1.64085319 #00FF00
    120 -0.18171090 #BEBEBE
    321  0.52274973 #DF5F5F
    324 -1.41307014 #00FF00
    325 -0.44015897 #5FDF5F
    326  0.30939426 #CE8F8F
    327  0.48981968 #DF5F5F
    328 -0.05552603 #BEBEBE
    329 -2.66488475 #00FF00
    330 -2.42844114 #00FF00
    331 -0.32863339 #8FCE8F
    332  0.68913235 #EF3030
    335 -0.48601589 #5FDF5F
    336 -0.01267029 #BEBEBE
    337 -0.23956608 #8FCE8F

    [[4]]$plot.data.cpd
    NULL


    [[5]]
    [[5]]$plot.data.gene
        kegg.names  labels
    18        4342     MOS
    19         983    CDK1
    20        4085  MAD2L1
    21        9133   CCNB2
    22        9088  PKMYT1
    23        6195 RPS6KA1
    24        5594   MAPK1
    25        5604  MAP2K1
    39         996   CDC27
    40         991   CDC20
    42        9232   PTTG1
    44         891   CCNB1
    45         983    CDK1
    46       26271   FBXO5
    47        6500    SKP1
    48         898   CCNE1
    49      286151  FBXO43
    50         996   CDC27
    51         991   CDC20
    52        9232   PTTG1
    53        6500    SKP1
    55         995  CDC25C
    56      286151  FBXO43
    57        1017    CDK2
    58        8945    BTRC
    59        8945    BTRC
    60         815  CAMK2A
    61        5347    PLK1
    62        6790   AURKA
    96         996   CDC27
    97         991   CDC20
    103      89869   PLCZ1
    105       3708   ITPR1
    107        801   CALM1
    113       5530  PPP3CA
    118      22849   CPEB3
    120      22849   CPEB3
    123       5347    PLK1
    129       6500    SKP1
    130       8945    BTRC
    131       5566  PRKACA
    132        107   ADCY1
    143       7529   YWHAB
    146       5499  PPP1CA
    148       1432  MAPK14
    150       5347    PLK1
    153       9700   ESPL1
    154       9700   ESPL1
    159       8243   SMC1A
    160       9985    REC8
    161      10734   STAG3
    164       8243   SMC1A
    165       9985    REC8
    166      10734   STAG3
    169     151648    SGO1
    170       5525 PPP2R5A
    171       5515  PPP2CA
    178       9126    SMC3
    179       5515  PPP2CA
    180       5241     PGR
    181        367      AR
    182       3480   IGF1R
    183       3630     INS
    184       3479    IGF1
    191     245711   SPDYA
    192        983    CDK1
    194       1017    CDK2
    198       9748     SLK
    207        699    BUB1
    208        995  CDC25C
    214       9126    SMC3
                                                                               all.mapped
    18                                                                                   
    19                                                                                983
    20                                                                    4085,8379,10459
    21                                                                               9133
    22                                                                               9088
    23                                                               6195,6196,6197,27330
    24                                                                          5594,5595
    25                                                                               5604
    39  996,8697,8881,10393,25847,25906,29882,29945,51433,51434,51529,64682,119504,246184
    40                                                                                991
    42                                                                               9232
    44                                                                                891
    45                                                                                983
    46                                                                              26271
    47                                                                     6500,8454,9978
    48                                                                           898,9134
    49                                                                             286151
    50  996,8697,8881,10393,25847,25906,29882,29945,51433,51434,51529,64682,119504,246184
    51                                                                                991
    52                                                                               9232
    53                                                                     6500,8454,9978
    55                                                                                995
    56                                                                             286151
    57                                                                               1017
    58                                                                         8945,23291
    59                                                                         8945,23291
    60                                                                    815,816,817,818
    61                                                                               5347
    62                                                                               6790
    96  996,8697,8881,10393,25847,25906,29882,29945,51433,51434,51529,64682,119504,246184
    97                                                                                991
    103                                                                                  
    105                                                                    3708,3709,3710
    107                                                          801,805,808,91860,163688
    113                                                          5530,5532,5533,5534,5535
    118                                                          22849,64506,80315,132864
    120                                                          22849,64506,80315,132864
    123                                                                              5347
    129                                                                    6500,8454,9978
    130                                                                        8945,23291
    131                                                                         5566,5567
    132                                                107,109,111,112,113,114,115,196883
    143                                                    7529,7531,7532,7533,7534,10971
    146                                                                    5499,5500,5501
    148                                                               1432,5600,5603,6300
    150                                                                              5347
    153                                                                              9700
    154                                                                              9700
    159                                                                              8243
    160                                                                              9985
    161                                                                             10734
    164                                                                              8243
    165                                                                              9985
    166                                                                             10734
    169                                                                            151648
    170                                                          5525,5526,5527,5528,5529
    171                                                               5515,5516,5518,5519
    178                                                                              9126
    179                                                               5515,5516,5518,5519
    180                                                                              5241
    181                                                                               367
    182                                                                              3480
    183                                                                                  
    184                                                                              3479
    191              245711,285955,387778,441272,441273,442590,729597,100310812,102723555
    192                                                                               983
    194                                                                              1017
    198                                                                              9748
    207                                                                               699
    208                                                                               995
    214                                                                              9126
        type    x   y width height    mol.data mol.col
    18  gene  258 319    46     17          NA #FFFFFF
    19  gene  437 713    46     17 -2.72622712 #00FF00
    20  gene  925 544    46     17 -1.94662076 #00FF00
    21  gene  437 697    46     17 -3.38030727 #00FF00
    22  gene  392 557    46     17 -2.70357324 #00FF00
    23  gene  392 362    46     17  2.50282218 #FF0000
    24  gene  392 319    46     17 -0.41712479 #5FDF5F
    25  gene  325 319    46     17 -0.25544203 #8FCE8F
    39  gene  653 610    46     17 -3.11519198 #00FF00
    40  gene  653 627    46     17 -3.23340373 #00FF00
    42  gene  653 683    46     17 -2.76607223 #00FF00
    44  gene  786 695    46     17 -2.66975910 #00FF00
    45  gene  786 712    46     17 -2.72622712 #00FF00
    46  gene  210 555    46     17 -1.41354365 #00FF00
    47  gene 1207 496    46     17 -0.73283181 #30EF30
    48  gene 1007 540    46     17 -1.61936153 #00FF00
    49  gene 1129 555    46     17 -1.90318480 #00FF00
    50  gene 1007 613    46     17 -3.11519198 #00FF00
    51  gene 1007 630    46     17 -3.23340373 #00FF00
    52  gene 1007 683    46     17 -2.76607223 #00FF00
    53  gene  210 488    46     17 -0.73283181 #30EF30
    55  gene  494 555    46     17 -2.88735131 #00FF00
    56  gene  703 555    46     17 -1.90318480 #00FF00
    57  gene 1007 557    46     17  0.83464695 #FF0000
    58  gene  210 505    46     17  1.32483774 #FF0000
    59  gene 1207 513    46     17  1.32483774 #FF0000
    60  gene 1129 466    46     17 -1.30920618 #00FF00
    61  gene 1178 466    46     17 -3.32923742 #00FF00
    62  gene  324 183    46     17 -2.02575149 #00FF00
    96  gene  255 610    46     17 -3.11519198 #00FF00
    97  gene  255 627    46     17 -3.23340373 #00FF00
    103 gene 1129 263    46     17          NA #FFFFFF
    105 gene 1130 341    46     17 -1.13400430 #00FF00
    107 gene 1129 418    46     17 -0.38777637 #8FCE8F
    113 gene 1066 466    46     17  1.16784574 #FF0000
    118 gene  278 237    46     17  0.48214326 #DF5F5F
    120 gene  786 295    46     17  0.48214326 #DF5F5F
    123 gene  842 249    46     17 -3.32923742 #00FF00
    129 gene  703 488    46     17 -0.73283181 #30EF30
    130 gene  703 505    46     17  1.32483774 #FF0000
    131 gene  157 257    46     17 -0.33879325 #8FCE8F
    132 gene  157 191    46     17  4.18255612 #FF0000
    143 gene  434 578    46     17 -2.66575328 #00FF00
    146 gene  561 483    46     17 -0.26271510 #8FCE8F
    148 gene  526 459    46     17  1.44870311 #FF0000
    150 gene  578 516    46     17 -3.32923742 #00FF00
    153 gene  653 722    46     17 -2.74549200 #00FF00
    154 gene 1007 722    46     17 -2.74549200 #00FF00
    159 gene  630 803    46     17 -0.64663966 #30EF30
    160 gene  653 769    46     17  0.54898951 #DF5F5F
    161 gene  653 786    46     17  0.17171287 #BEBEBE
    164 gene  984 803    46     17 -0.64663966 #30EF30
    165 gene 1007 769    46     17  0.54898951 #DF5F5F
    166 gene 1007 786    46     17  0.17171287 #BEBEBE
    169 gene  950 683    46     17 -2.65416038 #00FF00
    170 gene  950 722    46     17 -0.58635603 #5FDF5F
    171 gene  950 739    46     17 -2.19281051 #00FF00
    178 gene  676 803    46     17 -0.90701064 #00FF00
    179 gene 1094 490    46     17 -2.19281051 #00FF00
    180 gene  206 154    46     17  2.10326793 #FF0000
    181 gene  255 154    46     17 -0.28712304 #8FCE8F
    182 gene   98 154    46     17 -0.37082634 #8FCE8F
    183 gene   73 103    46     17          NA #FFFFFF
    184 gene  123 103    46     17  0.16957083 #BEBEBE
    191 gene  431 433    46     17  3.72750174 #FF0000
    192 gene  431 450    46     17 -2.72622712 #00FF00
    194 gene  487 436    46     17  0.83464695 #FF0000
    198 gene  578 556    46     17  0.07786989 #BEBEBE
    207 gene  860 500    46     17 -3.03236821 #00FF00
    208 gene 1204 636    46     17 -2.88735131 #00FF00
    214 gene 1030 803    46     17 -0.90701064 #00FF00

    [[5]]$plot.data.cpd
        kegg.names labels all.mapped     type    x   y width height mol.data
    26      C00410 C00410            compound  206 102     8      8       NA
    104     C01245 C01245            compound 1129 302     8      8       NA
    106     C00076 C00076            compound 1129 379     8      8       NA
    134     C00575 C00575            compound  157 225     8      8       NA
        mol.col
    26  #FFFFFF
    104 #FFFFFF
    106 #FFFFFF
    134 #FFFFFF

![](hsa03013.hsa03013.png) ![](hsa03030.hsa03030.png)
![](hsa03440.hsa03440.png)

![](hsa04110.hsa04110.png) ![](hsa04114.hsa04114.png)

\##Section 3. Gene Ontology (GO) We can also do a similar procedure with
gene ontology. Similar to above, go.sets.hs has all GO terms. go.subs.hs
is a named list containing indexes for the BP, CC, and MF ontologies.
Let’s focus on BP (a.k.a Biological Process) here.

``` r
library(gage)
library(gageData)

data(go.sets.hs)
data(go.subs.hs)

# Check the BP subset exists
names(go.subs.hs)  # should include "BP"
```

    [1] "BP" "CC" "MF"

``` r
gobpsets <- go.sets.hs[go.subs.hs$BP]

# Check foldchanges
head(foldchanges)           # must be named numeric vector
```

         1266     54855      1465      2034      2150      6659 
    -2.422719  3.201955 -2.313738 -1.888019  3.344508  2.392288 

``` r
all(names(foldchanges) %in% unlist(gobpsets))  # check overlap
```

    [1] FALSE

``` r
# Run GAGE
gobpres <- gage(foldchanges, gsets=gobpsets, same.dir=TRUE)

# Inspect results
lapply(gobpres, head)
```

    $greater
                                                 p.geomean stat.mean        p.val
    GO:0007156 homophilic cell adhesion       8.519724e-05  3.824205 8.519724e-05
    GO:0002009 morphogenesis of an epithelium 1.396681e-04  3.653886 1.396681e-04
    GO:0048729 tissue morphogenesis           1.432451e-04  3.643242 1.432451e-04
    GO:0007610 behavior                       1.925222e-04  3.565432 1.925222e-04
    GO:0060562 epithelial tube morphogenesis  5.932837e-04  3.261376 5.932837e-04
    GO:0035295 tube development               5.953254e-04  3.253665 5.953254e-04
                                                  q.val set.size         exp1
    GO:0007156 homophilic cell adhesion       0.1951953      113 8.519724e-05
    GO:0002009 morphogenesis of an epithelium 0.1951953      339 1.396681e-04
    GO:0048729 tissue morphogenesis           0.1951953      424 1.432451e-04
    GO:0007610 behavior                       0.1967577      426 1.925222e-04
    GO:0060562 epithelial tube morphogenesis  0.3565320      257 5.932837e-04
    GO:0035295 tube development               0.3565320      391 5.953254e-04

    $less
                                                p.geomean stat.mean        p.val
    GO:0048285 organelle fission             1.536227e-15 -8.063910 1.536227e-15
    GO:0000280 nuclear division              4.286961e-15 -7.939217 4.286961e-15
    GO:0007067 mitosis                       4.286961e-15 -7.939217 4.286961e-15
    GO:0000087 M phase of mitotic cell cycle 1.169934e-14 -7.797496 1.169934e-14
    GO:0007059 chromosome segregation        2.028624e-11 -6.878340 2.028624e-11
    GO:0000236 mitotic prometaphase          1.729553e-10 -6.695966 1.729553e-10
                                                    q.val set.size         exp1
    GO:0048285 organelle fission             5.841698e-12      376 1.536227e-15
    GO:0000280 nuclear division              5.841698e-12      352 4.286961e-15
    GO:0007067 mitosis                       5.841698e-12      352 4.286961e-15
    GO:0000087 M phase of mitotic cell cycle 1.195672e-11      362 1.169934e-14
    GO:0007059 chromosome segregation        1.658603e-08      142 2.028624e-11
    GO:0000236 mitotic prometaphase          1.178402e-07       84 1.729553e-10

    $stats
                                              stat.mean     exp1
    GO:0007156 homophilic cell adhesion        3.824205 3.824205
    GO:0002009 morphogenesis of an epithelium  3.653886 3.653886
    GO:0048729 tissue morphogenesis            3.643242 3.643242
    GO:0007610 behavior                        3.565432 3.565432
    GO:0060562 epithelial tube morphogenesis   3.261376 3.261376
    GO:0035295 tube development                3.253665 3.253665

\##Section 4. Reactome Analysis

(download reactome package)

Let’s now conduct over-representation enrichment analysis and
pathway-topology analysis with Reactome using the previous list of
significant genes generated from our differential expression results
above.

First, Using R, output the list of significant genes at the 0.05 level
as a plain text file:

``` r
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), "symbol"]
print(paste("Total number of significant genes:", length(sig_genes)))
```

    [1] "Total number of significant genes: 8147"

``` r
write.table(sig_genes, file="significant_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Q8: What pathway has the most significant “Entities p-value”? Do the
most significant pathways listed match your previous KEGG results? What
factors could cause differences between the two methods?

A8: The CellCycle has the most significant pathway, yes many major
pathways overlap like te ribosome, apoptosis ad metabolic pathways,
however reactome is more detailed and has more subdivided pathways

\##Section 5. GO online (OPTIONAL) Gene Set Gene Ontology (GO)
Enrichment is a method to determine over-represented or
under-represented GO terms for a given set of genes. GO terms are formal
structured controlled vocabularies (ontologies) for gene products in
terms of their biological function. The goal of this analysis is to
determine the biological process the given set of genes are associated
with.

Q9: What pathway has the most significant “Entities p-value”? Do the
most significant pathways listed match your previous KEGG results? What
factors could cause differences between the two methods?

A9:
