% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mGSZ2.R
\name{mGSZ2}
\alias{mGSZ2}
\title{mGSZ2: Improved Modified Gene Set Z-score method for gene set analysis.}
\usage{
mGSZ2(x, y, l, rankFn = "MA", min.sz = 5, pv = 0, w1 = 0.2,
  w2 = 0.5, vc = 10, p = 200, rankInParallel = !F)
}
\arguments{
\item{x}{Gene expression data matrix (rows as genes and columns as samples).
Raw counts for RNA-Seq data.}

\item{y}{Gene set data (dataframe/table/matrix/list)}

\item{l}{Vector of response values (example:1,2)}

\item{rankFn}{One of 'MA', 'RNA' if data comes from microarrays-compatible
or RNA-Seq technologies respectively. Or a function with inputs (x, l), and
returns a vector with the same length as the nrow(x).  Or a matrix of
rankings genes as rows (must have genes as rownames), first column for real
ranking, rest of columns for permuted rankings, in this case mGSZ2 will not
use x and l inputs (can be NA).}

\item{min.sz}{Minimum size of gene sets (number of genes in a gene set) to
be included in the analysis}

\item{pv}{Estimate of the variance associated with each observation}

\item{w1}{Weight 1, parameter used to calculate the prior variance
obtained with class size var.constant. This penalizes especially small
classes and small subsets. Default is 0.2. Values around 0.1 - 0.5
are expected to be reasonable}

\item{w2}{Weight 2, parameter used to calculate the prior variance
obtained with the same class size as that of the analyzed class. This
penalizes small subsets from the gene list. Default is 0.5. Values around
0.3 and 0.5 are expected to be reasonable}

\item{vc}{Size of the reference class used with wgt1. Default is 10}

\item{p}{Number of permutations for p-value calculation}

\item{rankInParallel}{If FALSE, permutation gene rankings will be calculated
sequentially. Useful if rankFn is provided and is already parallelized.}
}
\value{
Dataframe with gene sets (in decreasing order based on the
significance) reported by mGSZ method, scores, p-values, and list of genes
that contributed to the enrichment.
}
\description{
\code{mGSZ2} Gene set analysis based on Gene Set Z scoring function and
asymptotic p-value
}
\details{
A function for Gene set analysis based on Gene Set Z-scoring
function and asymptotic p- value. It differs from GSZ (Toronen et al 2009)
in that it implements asymptotic p-values instead of empirical p-values.
Asymptotic p-values are based on fitting suitable distribution model to
the permutation data. Unlike empirical p-values, the resolution of
asymptotic p-values are independent of the number of permutations and
hence requires consideralbly fewer permutations. In addition to GSZ, this
function allows the users to carry out analysis with seven other scoring
functions and compare the results.
}
\examples{
\dontrun{
data(dummyData);
mGSZres <- mGSZ2(dummyData$x, dummyData$y, dummyData$l);
}
}
\author{
Juan Cruz Rodriguez \email{jcrodriguez@bdmg.com.ar}
}
