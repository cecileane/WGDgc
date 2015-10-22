\name{sampleData1}
\alias{sampleData1}
\docType{data}

\title{Simulated gene count data with 1 WGD event}
\description{Sample gene count data simulated with 1 WGD, 4 species (A, B, C, D) and 6000 families.}

\usage{data(sampleData1)}

\format{
  A data frame with 6000 observations on the following 4 species as 4 named variables: A, B, C, D.
}

\details{
These data were generated according to the following species tree (in
simmap format version 1.1), with a single WGD event located on the
internal edge leading to the MRCA of species A and B and retention rate
0.6:
"(D:{0,18.03}, (C:{0,12.06},(B:{0,7.06},A:{0,7.06}):{0,2.50 :wgd,0:0,2.50}):{0, 5.97});"
The duplication and loss rates used for simulation were 0.02 and
0.03. Families with 0 or 1 copy were excluded. All families were started
with only one ancestral gene at the root of the species tree.
}

\examples{
data(sampleData1)
dat <- sampleData1[1:100,] # reducing data to run examples faster

tree1WGD.str = "(D:{0,18.03}, (C:{0,12.06},(B:{0,7.06},A:{0,7.06})
               :{0,2.50 :wgd,0:0,2.50}):{0, 5.97});"
# tree with a single WGD event along the edge to MRCA of species A and B 
tree1WGD = read.simmap(text=tree1WGD.str)

MLEGeneCount(tree1WGD, dat, dirac=1, conditioning="twoOrMore")
# to estimate retention, duplication and loss rates

MLEGeneCount(tree1WGD, dat, dirac=1, conditioning="twoOrMore",
             fixedRetentionRates=TRUE, startingQ=0.6)
# to estimate the duplication and loss rates only, 
# based on a hypothesized retention rate 0.6 at the WGD.

filtered <- subset(dat, (A>0| B>0 | C>0) & D>0 )
# families with at least one copy in both clades at the root

MLEGeneCount(tree1WGD, filtered,dirac=1,conditioning="oneInBothClades")
# uses the appropriate filtering

## Analysis under a tree with no WGD
tree0WGD.str = "(D:{0,18.03}, (C:{0,12.06},(B:{0,7.06},A:{0,7.06})
             :{0,5.00}):{0, 5.97});"
tree0WGD = read.simmap(text=tree0WGD.str)
MLEGeneCount(tree0WGD, dat, dirac=1, conditioning="twoOrMore", 
             fixedRetentionRates=TRUE)

## Analysis under a tree with 2 events: one WGD and one WGT
tree2events.str = "(D:{0,18.03}, (C:{0,12.06},(B:{0,7.06},A:{0,7.06}):
            {0,2.50 :wgt,0:0,2.50}):{0, 2.985: wgd,0:0,2.985});"
# oldest event: WGD on edge to MRCA of species A, B and C.
# recent event: WGT on edge to MRCA of species A, B
tree2events = read.simmap(text=tree2events.str)
MLEGeneCount(tree2events, dat, dirac=1, conditioning="twoOrMore")
}
