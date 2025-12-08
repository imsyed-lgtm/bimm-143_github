# Class 06: R Studio Functions
Iman Syed (A18596789)

All functions in R have at least 3 things:

- a **name**, we pick this and use it to call the function. \_ Input
  **arguements** , there can be multiple comma seperated inputs to the
  function.
- the **body** , lines of R code that do the work of the function.

Our first wee function

``` r
add <- function(x, y=1) {
  x+y
}
```

Let’s test our function

``` r
add( c(1,2,3), y=10)
```

    [1] 11 12 13

``` r
add(10,10)
```

    [1] 20

\##second function Let’s try something more interesting. Make a sequence
generation tool.

The ‘sample()’ function could be useful here.

``` r
sample(1:10, size= 3)
```

    [1] 8 6 9

Change this to work with nucleotides A C G and T and return 3 of them

``` r
n <- c("A", "C","G", "T")
sample(n,size=15, replace = TRUE)
```

     [1] "G" "T" "C" "G" "C" "A" "A" "C" "A" "G" "A" "G" "C" "C" "G"

Turn this snippet into a function that returns a user specified length
dna sequence. Let’s call it ‘generate_dna()’ …

``` r
generate_dna <- function(length=10, fasta= FALSE){
  n <- c("A", "C","G", "T")
  v <-sample(n,size=length, replace = TRUE)
  #Make a single element vector
  
  s <- paste(v, collapse="")
  cat("Well Done you! \n")
  
  if(fasta){
    return(s)
  } else {
    return (v)
  }

  generate_dna(5)
  }
```

``` r
s <- generate_dna(5)
```

    Well Done you! 

``` r
s
```

    [1] "A" "A" "A" "T" "C"

I want the option to return a single element character vector with my
sequence all together like this: “GGAGTAC”

``` r
generate_dna(10, fasta= FALSE)
```

    Well Done you! 

     [1] "G" "T" "A" "A" "C" "A" "C" "A" "A" "C"

Make a third function that generates a protein sequence of a user
unspecified length and format.

``` r
generate_protein <- function(size= 15, fasta = TRUE) {
  aa <- c(
    "A","R","N","D","C","E","Q","G","H","I",
    "L","K","M","F","P","S","T","W","Y","V")
  seq <- sample(aa, size = size, replace = TRUE)
  

  
  if (fasta) {
    return(paste(seq, collapse= ""))
  } else {
    return(seq)
  }
}
```

``` r
generate_protein(10)
```

    [1] "QRTAWYFSFG"

Q:Generate random protein sequences between lengths 5 and 12
amino-acids.

- One approach is to do this with brute force calling our function for
  each length between 5 and 12

``` r
generate_protein(5)
```

    [1] "TLKDT"

``` r
generate_protein(6)
```

    [1] "ARAKAY"

Another approach is to write a ‘for()’ loop to itterate over the input
valued 5 to 12.

``` r
seq_lengths <- 6:12
for (i in seq_lengths) {
  cat(">",i, "\n")
  cat(generate_protein(i))
  cat("\n")
  
}
```

    > 6 
    DHHHFL
    > 7 
    NCVWYLH
    > 8 
    TSHFFNIA
    > 9 
    VNIMDFNVP
    > 10 
    TAMNRNTVSA
    > 11 
    KMVQSANTNVE
    > 12 
    IWQKWTMNTTNK

-A very useful thrid R specific approach is to use the

``` r
sapply(6:12, generate_protein)
```

    [1] "TCNKKG"       "LLPCGGR"      "HRECAIIC"     "MYNMLSRPW"    "GMMYPNMAPD"  
    [6] "GCNEHRAKQCQ"  "VNHMMLHFDAIN"

**Key Point**: Writing functions in R is doable but not the easiest
thing. Starting with a working snippet of code and then using LLM tools
to improve and generalize your function code is a productive approach.
