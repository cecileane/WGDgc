\encoding{utf-8}
\name{logLik_CsurosMiklos}
\alias{logLik_CsurosMiklos}

\title{Log-likelihood of count data on a phylogenetic tree}

\description{Calculates the probability of gene count data on a phylogenetic tree under a birth-and-death process and whole genome duplication (or triplication) events, conditional on n surviving gene lineages at the root. Also computes the probability of a family going extinct.}

\usage{
logLik_CsurosMiklos(logLamlogMu, nLeaf, nFamily, phyloMat,
                    geneCountData, mMax, wgdTab, edgeOrder)
}

\arguments{
  \item{logLamlogMu}{vector of size 1 or 2, for the log of the
  duplication and loss rates. When a single rate is provided, the
  duplication and loss rates are assumed to be equal.}
  \item{nLeaf}{number of present-day species.}
  \item{nFamily}{number of gene families.}
  \item{phyloMat}{a phylogenetic matrix with 4 columns: parent (ancestor
  node), child (descendant node), time (branch length), and species
  names. The number of rows is the number of nodes in the tree.}
  \item{geneCountData}{data frame with one column for each species and
  one row for each family, containing the number of gene copies in each
  species for each gene family. The column names must match the species
  names in the tree.}
  \item{mMax}{maximum number of surviving lineages at the root, at which
  the likelihood will be computed.}
  \item{wgdTab}{data frame with 5 columns: node before event, event
  type (WGD or WGT) and retention rates of 1, 2 and 3 gene copies. The
  number of rows is the number of WGD events.}
  \item{edgeOrder}{a data frame listing the tree edges in post-order
  traversal with information on which are birth-death and WGD/T edges.}
}

\value{
 \item{loglikRoot}{matrix of size \code{nMax+1} by \code{nFamily} giving
 the log likelihood of each gene family given that there are n surviving
 gene lineages at the root in row n+1. Column k corresponds to family
 k.}
 \item{doomedRoot}{probability that a single gene lineage present at the
 root goes extinct.}
 \item{doomedRootLeft}{probability that a single gene lineage at the
 root goes extinct in the clade on the left side of the root.}
 \item{doomedRootRight}{probability that a single gene lineage at the
 root goes extinct in the clade on the right side of the root.}
}

\references{
Csuros M and Miklos I (2009). Streamlining and large ancestral genomes
in archaea inferred with a phylogenetic birth-and-death model.
\emph{Molecular Biology and Evolution}. \bold{26}:2087-2095.

Charles-Elie Rabier, Tram Ta and Cécile Ané (2013). Detecting and
locating whole genome duplications on a phylogeny: a probabilistic
approach. \emph{Molecular Biology and Evolution}. \bold{31}(3):750-762.
}

\author{Cécile Ané}

\seealso{
 \code{\link{processInput}, \link{getEdgeOrder}}.
}

\examples{
tre.string = "(D:{0,18.03},(C:{0,12.06},(B:{0,7.06},
              A:{0,7.06}):{0,2.49:wgd,0:0,2.50}):{0, 5.97});"
tre.phylo4d = read.simmap(text=tre.string)
dat = data.frame(A=c(2,2,3,1), B=c(3,0,2,1), C=c(1,0,2,2), D=c(2,1,1,1));
a = processInput(tre.phylo4d, startingQ=0.9)
logLik_CsurosMiklos(log(c(.01,.02)), nLeaf=4, nFamily=4,
                    a$phyloMat,dat,mMax=8,a$wgdTab, a$edgeOrder)
}
