# Lab 5: Data Viz with ggplot
Iman Syed (A18596789)

Today we are exploring the “ggplot” package and how to make nice figures
in R.

There are lots of ways to make figures and plots in R.These include:-
-so called “base” R - and add a package like **ggplot2**

here is a simple “base” R plot.

``` r
options(ggrepel.max.overlaps = Inf)

head(cars)
```

      speed dist
    1     4    2
    2     4   10
    3     7    4
    4     7   22
    5     8   16
    6     9   10

We can simply pass to the \`plot()~ function

``` r
plot(cars)
```

![](lab05_files/figure-commonmark/unnamed-chunk-2-1.png)

Key point: Base R is quick but not so nice and simple looking in some
folks eyes

Lets see how we can plot this with **ggplot2**

``` r
library(ggplot2)
ggplot(cars)
```

![](lab05_files/figure-commonmark/unnamed-chunk-3-1.png)

So an error may pop up, so Ist you need to install an add-on package.
for this we use the `install.packages()` function. however, we do this
in the CONSOLE, NOT our report\*\*

2nd we need to load the package with the `library()` function everytime
we want to use it.

``` r
ggplot(cars)
```

![](lab05_files/figure-commonmark/unnamed-chunk-4-1.png)

Every ggplot is composed of at least 3 layers:

\-**data** (i.e a data.frame with the things you want to plot),
-aesthetics **aes()** that map the columns of data to your plot features
-geoms like **geom_point()** that sort hoe the plot appears

``` r
ggplot(cars)+
  aes(x= speed, y=dist) +
  geom_abline()
```

![](lab05_files/figure-commonmark/unnamed-chunk-5-1.png)

Lest add more layers to our ggplot add a line shoing the relationship
between x and y

``` r
ggplot(cars)+
  aes(x=speed, y=dist)+
  geom_point()+   #this causes points to show on the graph
  geom_smooth (method ="lm", se= FALSE) +   # geom_smooth adds a trendline, (method ="lm", se= FALSE) makes it straight
  labs(title="Silly plot of Speed vs Stopping distance")
```

    `geom_smooth()` using formula = 'y ~ x'

![](lab05_files/figure-commonmark/unnamed-chunk-6-1.png)

Going Further:

Read some genomic data into the report

``` r
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

Q1:How many genes are in the data set?

``` r
nrow(genes)   #counts how many rows in a data set
```

    [1] 5196

``` r
ncol(genes)   #counts how many rows are in a data set
```

    [1] 4

Q2: how many up regualated genes are there?

``` r
sum(genes$State == "up")
```

    [1] 127

a useful function for counting up occurrances of things in a vector is
the table function

``` r
table (genes$State)
```


          down unchanging         up 
            72       4997        127 

Make a v1 figure:

# Template for plots:

``` r
p <- ggplot(genes) +
  aes(x= Condition1,y=Condition2, col= State)+   # the col=State function gives your plots colours
  geom_point()
  
  
p
```

![](lab05_files/figure-commonmark/unnamed-chunk-12-1.png)

Now we can customize the colors using this function:

``` r
p + scale_colour_manual( values=c("purple","green","cyan") )
```

![](lab05_files/figure-commonmark/unnamed-chunk-13-1.png)

After changing the colors, we can add labels and more layers:

``` r
p + scale_colour_manual(values=c("purple","green","cyan")) +
    labs(title="Gene Expresion Changes Upon Drug Treatment",
         x="Control (no drug) ",
         y="Drug Treatment")
```

![](lab05_files/figure-commonmark/unnamed-chunk-14-1.png)

Section 7: More Plotting examples

Read in gapminder dataset:

``` r
# File location online
url <- "https://raw.githubusercontent.com/jennybc/gapminder/master/inst/extdata/gapminder.tsv"

gapminder <- read.delim(url)
```

Lets take a little look at this dataset

``` r
head(gapminder, 3)  #top of the data set in alphabetical order
```

          country continent year lifeExp      pop gdpPercap
    1 Afghanistan      Asia 1952  28.801  8425333  779.4453
    2 Afghanistan      Asia 1957  30.332  9240934  820.8530
    3 Afghanistan      Asia 1962  31.997 10267083  853.1007

``` r
tail( gapminder, 3) #tail of the data set in alphabetical order
```

          country continent year lifeExp      pop gdpPercap
    1702 Zimbabwe    Africa 1997  46.809 11404948  792.4500
    1703 Zimbabwe    Africa 2002  39.989 11926563  672.0386
    1704 Zimbabwe    Africa 2007  43.487 12311143  469.7093

Q4: How many different country values are in this dataset?

``` r
# nrow() does not work here, because many countries are repeated in rows!
length(table(gapminder$country))
```

    [1] 142

Q5:How many different continents values are in this dataset?

``` r
#we can use the same way as Q4, however you can use a more efficient unique() function

unique(gapminder$continent)
```

    [1] "Asia"     "Europe"   "Africa"   "Americas" "Oceania" 

Making a rough plot of our data

``` r
ggplot(gapminder)+
  aes(gdpPercap,lifeExp)+
  geom_point()
```

![](lab05_files/figure-commonmark/unnamed-chunk-20-1.png)

Making the graphics prettier!!

``` r
ggplot(gapminder)+
  aes(gdpPercap,lifeExp,
  col=continent, label=country,)+
  geom_text()
```

![](lab05_files/figure-commonmark/unnamed-chunk-21-1.png)

I can use the **ggrepel** package to make more sensible lables here.

``` r
library(ggrepel)

ggplot(gapminder)+
  aes(gdpPercap, lifeExp, col= continent, label= country) +
  geom_point()+ 
  geom_text_repel() # only loads text that is not overlapping
```

![](lab05_files/figure-commonmark/unnamed-chunk-22-1.png)

I want a separate pannel per continent:

``` r
library(ggrepel)

ggplot(gapminder)+
  aes(gdpPercap, lifeExp, col= continent, label= country) +
  geom_point()+
  facet_wrap(~continent)
```

![](lab05_files/figure-commonmark/unnamed-chunk-23-1.png)

Q6: What are the main advantages of ggplot over base R plot?
