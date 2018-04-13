# mGSZ2
### Improved Modified Gene Set Z-score method for gene set analysis.
mGSZ2 is an improved version of the functional analysis [mGSZ R library](https://cran.r-project.org/web/packages/mGSZ/).

### Improvements
* Runs multicore by implementing the [BiocParallel library](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html).
* Allows the use of any desired genes ranking function.
* As output it returns the gene set score. Positive or negative depending on enriching condition.
* As output it returns the genes list that contributed to enrich each gene set, i.e., genes conforming the enrichment score leading edge.

### Installation
In R console type:

    R > install.packages('devtools');
    R > library(devtools);
    R > install_github('jcrodriguez1989/mGSZ2');

### Usage
    R > library('mGSZ2');
    R >
    R > # Load example data
    R > data(dummyData); # list that contains 'x', 'y', 'l'.
    R > # x is the expression matrix (genes as rows, samples as cols), y are the gene sets
    R > # l are the conditions (one for each sample from x).
    R >
    R > set.seed(8818);
    R > mGSZres <- mGSZ2(dummyData$x, dummyData$y, dummyData$l);
    R > head(mGSZres);
    gene.sets     pvalue mGszScore                                        impGenes
    set9       set9 0.05117402  1.952384 g55, g70, g68, g30, g51, g81, g44, g88, g75, g5
    set8       set8 0.06463807 -2.146767     g50, g92, g71, g54, g10, g46, g87, g13, g48
    set16     set16 0.10781461  1.907799 g53, g22, g68, g2, g18, g32, g80, g15, g63, g25
    set2       set2 0.12166137  2.124259                              g1, g61, g100, g67
    set1       set1 0.13451610 -1.879446                                    g4, g50, g11
    set7       set7 0.15552483  1.503978                                         g1, g76
