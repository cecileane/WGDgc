\encoding{utf-8}
\name{rgenecount}
\alias{rgenecount}

\title{Random generation of family sizes}

\description{Generates gene count data for multiple families
  along a phylogeny, using background rates of duplication and loss
  and possible whole genome duplication (WGD) or triplication (WGT)
  event(s), each with its own retention rate.}

\usage{
rgenecount(nfam, tre, lambdamu, retention, geomMean=NULL, dirac=NULL, 
           conditioning=c("none"))
}

\arguments{
  \item{nfam}{number of families to simulate}
  \item{tre}{a species tree in SIMMAP format.}
  \item{lambdamu}{vector of size 1 or 2, for the duplication rate
    (\eqn{\lambda}{lambda}) and loss rate (\eqn{\mu}{mu}).
    A vector of size 1 sets \eqn{\lambda=\mu}{lambda=mu}.}
  \item{retention}{vector of length the number of WGD/WGT events in the
  tree, giving the retention rate at each event. For WGTs, both extra
  copies are assumed to be retained independently with the same retention rate.}
  \item{geomMean}{the mean of the prior geometric distribution for the
  number of genes at the root.}
  \item{dirac}{value for the number of genes at the root, if
  fixed to the same value for all families.}
  \item{conditioning}{type of filtering. No filterning implemented yet.} 
}

\value{
  matrix with \code{nfam} rows, one per simulated family,
  and one column per node in the tree (tips and internal nodes).
}

\details{
  For the simmap format, see \code{\link{MLEGeneCount}}
  For WGT events, the 2 extra copies are assumed to be retained
  independently. With retention rate \eqn{q}, the probability to retain
  all 3 gene copies is then \eqn{q^2}, the probability to retain 2 gene
  copies is \eqn{2*q*(1-q)}, and the probability to retain the original
  gene only is \eqn{(1-q)^2}.

  The \code{geomMean} and \code{dirac} options are incompatible. 
}

\author{Cécile Ané}

\examples{
# tree with 2 WGDs. The second is placed immediately after
#                   the split between C and AB:
tre.string <- "(D:{0,18.03},(C:{0,12.06},(B:{0,7.06},
  A:{0,7.06}):{0,2.49:wgd,0:0,2.50:wgd,0:0,1e-10}):{0, 5.97});"
tre.phylo4d = read.simmap(text=tre.string)
# do this to see how edges and nodes are numbered,
#                which WGD is the first, which is the second:
processInput(tre.phylo4d, startingQ=c(.6,.2))
rgenecount(nfam=10,tre.phylo4d,lambdamu=c(.03,.04),retention=c(.6,.2),dirac=1)
}
