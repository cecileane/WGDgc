\name{sampleData2}
\alias{sampleData2}
\docType{data}

\title{Simulated gene count data with two WGD events}
\description{Sample gene count data simulated with 2 WGDs on the same branch, 4 species (A, B, C, D) and 6000 families.}

\usage{data(sampleData2)}

\format{
  A data frame with 6000 observations on the following 4 species as 4 named variables: A, B, C, D.
}

\details{
These data were generated according to the following species tree (in simmap format version 1.1), with both WGD events located along the internal edge leading species D, with retention rate 0.6 for the oldest event and 0.2 for the most recent event:

"(D:{0,6.01:0.2,0:0,6.01:0.6,0:0,6.01}, (C:0,12.06,(B:0,7.06,A:0,7.06):0,4.99):0,5.97);"

The duplication and loss rates used for simulation were 0.02 and 0.03. Families with 0 or 1 copy were excluded. All families were started with only one ancestral gene at the root of the species tree.
}

\examples{
data(sampleData2)
dat <- sampleData2[1:200,] # reducing data to run examples faster

tree2WGD.str="(D:{0,6.01:wgd,0:0,6.01:wgd,0:0,6.01}, (C:{0,12.06},
                (B:{0,7.06},A:{0,7.06}):{0,4.99}):{0,5.97});"
# both WGD events are located on the edge leading to species D
tree2WGD = read.simmap(text=tree2WGD.str)
MLEGeneCount(tree2WGD, dat, dirac=1, conditioning="twoOrMore")
}
