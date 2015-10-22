\encoding{utf-8}
\name{getLikGeneCount}
\alias{getLikGeneCount}

\title{Negative log-likelihood of gene count data}

\description{Calculates the overall negative log-likelihood of gene
  count data on a phylogenetic tree under a birth-and-death process and
  whole genome duplication events.}

\usage{
getLikGeneCount(para, input, geneCountData, mMax=NULL,
                geomProb=NULL, dirac=NULL, useRootStateMLE=FALSE,
                conditioning=c("oneOrMore", "twoOrMore",
                "oneInBothClades", "none"),
                equalBDrates=FALSE, fixedRetentionRates=TRUE)
}

\arguments{
  \item{para}{vector of parameters (see Details)}
  \item{input}{object output by function \code{processInput}}
  \item{geneCountData}{data frame with one column for each species and
  one row for each family, containing the number of gene copies in each
  species for each gene family. The column names must match the species
  names in the tree.}
  \item{mMax}{maximum number of surviving lineages at the root, at which
  the likelihood will be evaluated.}
  \item{geomProb}{inverse of the prior mean number of gene lineages at
  the root.}
  \item{dirac}{value for the number of genes at the root, when this is
  assumed to have a fixed value (according to a dirac prior
  distribution).}
  \item{useRootStateMLE}{if \code{TRUE}, the most likely number of genes
  at the root is determined for each family separately and is used to
  evaluate the likelihood function.}
  \item{conditioning}{type of conditioning for the likelihood
  calculation. The default is to calculate conditional probabilities on
  observing families with at least 1 gene copy (see Details in
  \code{\link{MLEGeneCount}}).}
  \item{equalBDrates}{if \code{TRUE}, the duplication and loss rates are equal.}
  \item{fixedRetentionRates}{if \code{TRUE}, it uses retention rates
  present in \code{input$wgdTab}. If \code{FALSE}, it uses retention
  rates in \code{para}.}
}

\details{
The vector \code{para} for the parameters to be
used is of size 1+number of WGD/Ts if the birth and death rates are
assumed equal, or 2+number of WGD/Ts otherwise. It starts with
\code{log(StartingBDrates[1])} if \code{equalBDrates} is \code{TRUE},
with \code{log(StartingBDrates)} otherwise. The remaining components
correspond to retention rates.
}

\value{
negative log-likelihood value
}

\references{
Csuros M and Miklos I (2009). Streamlining and large ancestral genomes
in archaea inferred with a phylogenetic birth-and-death model.
Molecular Biology and Evolution. \bold{26}:2087-2095.

Charles-Elie Rabier, Tram Ta and Cécile Ané (2013). Detecting and
locating whole genome duplications on a phylogeny: a probabilistic
approach. Molecular Biology and Evolution. \bold{31}(3):750-762.
}

\seealso{
 \code{\link{MLEGeneCount}, \link{logLik_CsurosMiklos}}.
}

\examples{
tre.string = "(D:{0,18.03},(C:{0,12.06},(B:{0,7.06},
              A:{0,7.06}):{0,2.49:wgd,0:0,2.50}):{0, 5.97});"
tre.phylo4d = read.simmap(text=tre.string)
dat = data.frame(A=c(2,2,3,1), B=c(3,0,2,1), C=c(1,0,2,2), D=c(2,1,1,1));
a = processInput(tre.phylo4d, startingQ=0.9)
getLikGeneCount(log(c(.01,.02)),a,dat,mMax=8,geomProb=1/1.5,
                conditioning="oneOrMore")
}

