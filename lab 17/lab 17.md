# lab 17
iman syed 18596789

##Section 1: Identify genetic variants of interest

There are a number of gene variants associated with childhood asthma, .
A study from Verlaan et al. (2009) shows that 4 candidate SNPs
demonstrate significant evidence for association. You want to find out
what they are by visiting OMIM (http://www.omim.org) and locating the
Verlaan et al. paper description.

Here we read the CSV file drived from Ensembl

``` r
mxl <- read.csv("MXL.csv")
head(mxl)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  NA19648 (F)                       A|A ALL, AMR, MXL      -
    2                  NA19649 (M)                       G|G ALL, AMR, MXL      -
    3                  NA19651 (F)                       A|A ALL, AMR, MXL      -
    4                  NA19652 (M)                       G|G ALL, AMR, MXL      -
    5                  NA19654 (F)                       G|G ALL, AMR, MXL      -
    6                  NA19655 (M)                       A|G ALL, AMR, MXL      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

``` r
table(mxl$Genotype..forward.strand.)/nrow(mxl)*100
```


        A|A     A|G     G|A     G|G 
    34.3750 32.8125 18.7500 14.0625 

Let us look at a separate population:

Read in the file for GBR

``` r
gbr <- read.csv("GBR.csv")
head(gbr)
```

      Sample..Male.Female.Unknown. Genotype..forward.strand. Population.s. Father
    1                  HG00096 (M)                       A|A ALL, EUR, GBR      -
    2                  HG00097 (F)                       G|A ALL, EUR, GBR      -
    3                  HG00099 (F)                       G|G ALL, EUR, GBR      -
    4                  HG00100 (F)                       A|A ALL, EUR, GBR      -
    5                  HG00101 (M)                       A|A ALL, EUR, GBR      -
    6                  HG00102 (F)                       A|A ALL, EUR, GBR      -
      Mother
    1      -
    2      -
    3      -
    4      -
    5      -
    6      -

Find proporton of G|G

``` r
round(table(gbr$Genotype..forward.strand.)/nrow(gbr)* 100, 2 )
```


      A|A   A|G   G|A   G|G 
    25.27 18.68 26.37 29.67 

The variant that is associated with childhood asthma is more frequent in
GBR populations than the MKL population.

##Section 4 :Population scale analysis

One sample is obviously not enough to know what is happening in a
population. You are interested in assessing genetic differences on a
population scale. So, you processed about ~230 samples and did the
normalization on a genome level. Now, you want to find whether there is
any association of the 4 asthma-associated SNPs (rs8067378…) on ORMDL3
expression.

Q13: Read this file into R and determine the sample size for each
genotype and their corresponding median expression levels for each of
these genotypes:

``` r
expres <-read.table("rs8067378_ENSG00000172057.6.txt")
head(expres)
```

       sample geno      exp
    1 HG00367  A/G 28.96038
    2 NA20768  A/G 20.24449
    3 HG00361  A/A 31.32628
    4 HG00135  A/A 34.11169
    5 NA18870  G/G 18.25141
    6 NA11993  A/A 32.89721

``` r
nrow(expres)
```

    [1] 462

A13: a) There are 462 samples

``` r
table(expres$geno)
```


    A/A A/G G/G 
    108 233 121 

A13: b) median expression levels for each of these genotype A/A 108 A/G
233 G/G 121

``` r
library(ggplot2)
```

let’s make a plot with colours:

Q14: Generate a boxplot with a box per genotype, what could you infer
from the relative expression value between A/A and G/G displayed in this
plot? Does the SNP effect the expression of ORMDL3?

``` r
ggplot(expres)+ aes(geno,exp,fill=geno)+
  geom_boxplot(notch= TRUE)
```

![](lab-17.markdown_strict_files/figure-markdown_strict/unnamed-chunk-9-1.png)

A14: We can see that this location the relative expression value of the
A/A genotype is the highest whilst G/G has the lowest expression value.
This suggests that this SNP does affect the expression of the ORMDL3.
