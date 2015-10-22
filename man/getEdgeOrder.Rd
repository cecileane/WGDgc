\encoding{utf-8}
\name{getEdgeOrder}
\alias{getEdgeOrder}

\title{list the tree nodes in a post-order traversal}

\description{Preprocessing to list the edges in a post-order traversal,
  for future use in likelihood calculation. The output includes
  information on which edges the birth-death process applies to, and
  which edges represent a whole genome duplication or triplication event.}

\usage{
getEdgeOrder(phyloMat,nLeaf,wgdTab)
}

\arguments{
  \item{phyloMat}{Matrix representation of the species tree and WGD events}
  \item{nLeaf}{Number of present-day species (i.e. number of leaves)}
  \item{wgdTab}{Table representation of WGD events with retention rates}
}

\value{
Data frame listing the edges in a post-order traversal, with the following components
  \item{child}{index of the edge's child node}
  \item{edge}{index of the edge, i.e. its row in phyloMat}
  \item{type}{"BD" if birth-death edge, "WGD" or "WGT" if the edge is modelling a WGD/T event, or "rootPrior" if the edge is parent to the root node}
  \item{scdsib}{\code{TRUE} if the edge is listed after a sibling edge, \code{FALSE} otherwise}
}

\details{
This function assumes that speciation nodes in phyloMat are given lower indices than singleton nodes when the tree is read in by phyext, that speciation nodes are in pre-order in \code{phyloMat}, and that 2 singleton nodes are used to represent each WGD.
}

\author{Cécile Ané}

\seealso{
 \code{\link{processInput}}.
}

