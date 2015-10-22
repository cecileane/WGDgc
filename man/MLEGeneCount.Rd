\encoding{utf-8}
\name{MLEGeneCount}
\alias{MLEGeneCount}

\title{Maximum likelihood estimation of gene turnover rates with WGD}

\description{Uses gene count data to estimates rates of gene duplication
  and gene loss along a phylogeny with zero, one or more whole genome
  duplication (WGD) or triplication (WGT) events. Also estimates the
  gene retention rate after each WGD/WGT event.}


\usage{

MLEGeneCount(tr, geneCountData, mMax=NULL, geomMean=NULL,
             dirac=NULL, useRootStateMLE=FALSE,
             conditioning=c("oneOrMore", "twoOrMore",
             "oneInBothClades", "none"),
             equalBDrates=FALSE, fixedRetentionRates=FALSE,
             startingBDrates=c(0.01, 0.02),startingQ=NULL)
}

\arguments{
  \item{tr}{a species tree in SIMMAP format (see Details).}
  \item{geneCountData}{data frame with one column for each species and
  one row for each family, containing the number of gene copies in each
  species for each gene family. The column names must match the species
  names in the tree.}
  \item{mMax}{maximum number of surviving lineages at the root, at which
  the likelihood will be computed.}
  \item{geomMean}{the mean of the prior geometric distribution for the
  number of genes at the root.}
  \item{dirac}{value for the number of genes at the root, when this is
  assumed to have a fixed value (according to a dirac prior
  distribution).}
  \item{useRootStateMLE}{if \code{TRUE}, the most likely number of
  surviving genes at the root is determined for each family separately,
  and is used to calculate the overall likelihood of the data. This
  value at the root may vary with the parameter values during likelihood
  optimization.}
  \item{conditioning}{type of conditioning for the likelihood
  calculation. The default is to calculate conditional probabilities on
  observing families with at least 1 gene copy (see Details).} 
  \item{equalBDrates}{if \code{TRUE}, the duplication and loss rates are
  constrained to be equal.}
  \item{fixedRetentionRates}{if \code{TRUE}, retention rates from the
  user-defined tree are fixed and used as provided. If \code{FALSE},
  retention rates are considered as parameters and are estimated by
  maximum likelihood.}
  \item{startingBDrates}{Vector of size 2, for the starting values of
  the duplication and loss rates. When \code{equalBDrates=TRUE}, only
  the first component is used.}
  \item{startingQ}{Vector of starting values for the retention rates at
  the WGD and WGT events.}
}

\value{
  \item{birthrate}{birth or duplication rate}
  \item{deathrate}{death or loss rate}
  \item{loglikelihood}{log of the likelihood}
  \item{WGDtable}{a WGD table with 5 columns: node before WGD/WGT, event
  type, and probabilities that 1, 2 or 3 gene copies are retained. The
  number of rows is the number of WGD/WGT events.}
  \item{phyloMat}{data frame with 5 columns to describe the phylogeny:
  parent (ancestor node), child (descendant node), time (branch length),
  species names and edge type (e.g. BD or WGD). The number of rows is
  the number of nodes in the tree.}
  \item{call}{initial call to the function}
  \item{convergence}{optimization convergence flag from the \code{optim}
  call. 0 means successful convergence.}
  \item{mMax}{mMax value used for the likelihood calculations}
}

\details{
The tree needs to be in simmap format (version 1.1). 
This format is similar to the newick parenthetical format, except that
branch lengths are given inside brackets where states are indicated at
specific times along each branch. Along a given branch, the token "0,18"
indicates state 0 for a duration of 18 time units. Tokens are separated
with ":". State 0 is used to indicate branch segments where only the
birth/death process applies for gene duplications and losses. Labels
"wgd" or "WGT" are used for branch segments at WGD events, and "wgt" or
"WGT" for segments at WGT events. Such segments need to have a length of
0.

For WGT events, the 2 extra copies are assumed to be retained
independently. With retention rate \eqn{q}, the probability to retain
all 3 gene copies is then \eqn{q^2}, the probability to retain 2 gene
copies is \eqn{2*q*(1-q)}, and the probability to retain the original
gene only is \eqn{(1-q)^2}.

Four types of conditional likelihoods are implemented. The option
\code{conditioning} should match the data filtering process: use
\code{conditioning="oneOrMore"} if all families with one or more gene
copies are included in the data, use \code{"twoOrMore"} to condition on
families having two of more genes, \code{"oneInBothClades"} if the data
set was filtered to include only families with at least one gene copy in
each of the two main clades stemming from the root. Unconditional
likelihoods are used with \code{conditioning="none"}.

The \code{geomMean}, \code{dirac} and \code{useRootStateMLE} options are incompatible. 

By default, \code{mMax} is set to the maximum family size for an exact
likelihood calculation. For data sets with one or more very large
families, this can cause \code{mMax} to be very large and calculation to
be very slow. In such cases, the user can set \code{mMax} to a lower
value to speed up calculations, at the cost of an approximation to the
likelihood of families with a larger family size.
}

\references{
Bailey, N. (1964) The Elements of Stochastic Processes. New York: John Wiley \& Sons

Bollback J. P. (2006) SIMMAP: Stochastic character mapping of discrete
traits on phylogenies. \emph{Bioinformatics}. \bold{7}:88 

De Bie, T. and Cristianini, N. and Demuth, J.P. and Hahn, M.W. (2006)
CAFE: a computational tool for the study of gene family
evolution. \emph{Bioinformatics}. \bold{22}:1269--1271

Hahn, M.W. and De Bie, T. and Stajich, J.E. and Nguyen, C. and
Cristianini, N. (2005) Estimating the tempo and mode of gene family
evolution from comparative genomic data.
\emph{Genome Res.}. \bold{15}:1153--1160

Crawford, F., Suchard, M. (2012) Transition probabilities for general
birth-death processes with applications in ecology, genetics, and evolution.
\emph{J Math Biol}. \bold{65}:553-580

Rabier, C., Ta, T. and Ané, C. (2013) Detecting and Locating Whole
Genome Duplications on a phylogeny: a probabilistic
approach. \emph{Molecular Biology and Evolution}. \bold{31}(3):750-762.
}

\author{Tram Ta, Charles-Elie Rabier}

\seealso{
  \code{\link{sampleData1}, \link{sampleData2}} for more examples.
}

\examples{
tre.string = "(D:{0,18.03},(C:{0,12.06},(B:{0,7.06},
              A:{0,7.06}):{0,2.49:wgd,0:0,2.50}):{0, 5.97});"
# tree with a single hypothesized WGD event, along the
# internal edge leading to the MRCA of species A and B
tre.phylo4d = read.simmap(text=tre.string)
tre.phylo   = as(tre.phylo4d, "phylo")
\dontrun{plot(collapse.singles(tre.phylo))}
dat = data.frame(A=c(2,2,3,1), B=c(3,0,2,1), C=c(1,0,2,2), D=c(2,1,1,1))
MLEGeneCount(tre.phylo4d, dat, geomMean=1.5, 
             conditioning="oneOrMore", fixedRetentionRates=TRUE)
}
