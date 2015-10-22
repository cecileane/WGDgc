\name{processInput}
\alias{processInput}

\title{Preprocessing function}

\description{Checking arguments and preparing data for future optimization}

\usage{
processInput(tr, equalBDrates=FALSE, fixedRetentionRates=TRUE,
             startingBDrates=c(0.01, 0.02),startingQ=NULL)
}

\arguments{
  \item{tr}{a species tree in SIMMAP format (see Details of function
  MLEGeneCount).}
  \item{equalBDrates}{if \code{TRUE}, the duplication and loss rates are
  equal.}
  \item{fixedRetentionRates}{if \code{TRUE}, retention rates will be
  fixed to \code{startingQ} during the future optimization. If
  \code{FALSE}, retention rates will be considered as parameters and
  will be estimated by maximum likelihood.}
  \item{startingBDrates}{Vector of size 2 as starting values for the
  duplication and loss rates. When \code{equalBDrates=TRUE} only the
  first component is used.}
  \item{startingQ}{Vector of starting values for retention
  rates. Default is 0.5 for all WGD events.}
}

\value{
  \item{phyloMat}{data frame to represent the phylogeny. The number of rows
  is the number of nodes in the species tree. There are 5 columns
  (Parent, Child, Time, Species, type).}
  \item{nLeaf}{number of present-day species (i.e. number of leaves)}
  \item{nNode}{number of nodes in the species tree}
  \item{wgdTab}{data frame with 5 columns. Each row
  corresponds to a WGD(s) or WGT(s). The first column
  gives the node just before the WGD/T. The second column 'type' says if
  the event is a WGD or WGT. The remaining columns contain the
  probabilities that only the original gene is retained, or if 2 (or 3)
  gene copies are retained.}
  \item{para}{Vector of parameters to be optimized. see Details}
  \item{lower}{Lower bounds for later optimization. see Details}
  \item{upper}{Upper bounds for later optimization. see Details}
}

\details{
The vector \code{para} of starting values for the parameters to be
optimized is of size 1+number of WGDs if the birth and death rates are
assumed equal, or 2+number of WGDs otherwise. It starts with
\code{log(StartingBDrates[1])} if \code{equalBDrates} is \code{TRUE},
with \code{log(StartingBDrates)} otherwise, and the remaining components
(corresponding to the retention rates) are \code{startingQ} if
\code{startingQ} is provided, 0.5 otherwise.

For WGT events, the 2 extra copies are assumed to be retained
independently. With retention rate \eqn{q}, the probability to retain
all 3 gene copies is then \eqn{q^2}, the probability to retain 2 gene
copies is \eqn{2*q*(1-q)}, and the probability to retain the original
gene only is \eqn{(1-q)^2}.

\code{lower} and \code{upper} are vectors whose sizes correspond to the
number of parameters for the lower and upper bounds of the different
parameters in a subsequent optimization search. The log of the
duplication and loss rates are unconstrained, while duplicate retention
rates are constrained in [0,1].
}

\examples{
tre.string = "(D:{0,18.03},(C:{0,12.06},(B:{0,7.06},
              A:{0,7.06}):{0,2.49:wgd,0:0,2.50}):{0, 5.97});"
tre.phylo4d = read.simmap(text=tre.string)
processInput(tre.phylo4d)
}
