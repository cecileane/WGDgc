MLEGeneCount <- function(tr, geneCountData, mMax=NULL, geomMean=NULL,
                         dirac=NULL, useRootStateMLE=FALSE,
                         conditioning=c("oneOrMore", "twoOrMore", "oneInBothClades", "none"),
                         equalBDrates=FALSE, fixedRetentionRates=FALSE,
                         startingBDrates=c(0.01, 0.02),startingQ=NULL){
  ## tr has simmap format, i.e. phylo4d objects from read.simmap()
  ## species names on geneCountData should match species names (leaves) on the tree
  if (isTRUE(all.equal(geneCountData, round(geneCountData))))
    geneCountData<-round(geneCountData) # check for integer data
  else
    stop("data are not of integer type")
  if (dim(geneCountData)[2]!=length(tipLabels(tr))) # check species	
    stop("the number of species in the tree does not match the number of species in the data file")
  if (sum(names(geneCountData)%in% tipLabels(tr)) != length(tipLabels(tr)))
    stop("the column names do not match the species names")
  if(is.null(mMax)) # initialize mMax to max # of surviving gene lineages
    mMax = max(apply(geneCountData,1,sum))
  if ((mMax!=floor(mMax)) || mMax<=0)
    stop("mMax must be a positive integer. Was ",mMax)
  if (mMax< max(apply(geneCountData,1,sum)))
    warning("the value of the parameter mMax should be larger to avoid approximations.",
            immediate.=TRUE)
  if (mMax< max(geneCountData))  
    stop(paste("the value of mMax should be no smaller than the largest gene count:",
               max(geneCountData)))

## old code with nPos: max possible number of genes at each internal node of the phylogeny.
# nPos = max(max(geneCountData)*3 +1, 101)
# if (nPos<(max(geneCountData)+1)|nPos<0)
#  stop("nPos is too small. It should be at least three times 
# the largest value observed in present-day species + 1. ")
# if (nPos<(max(geneCountData)*3+1))
#  warning("the value of the parameter nPos should be larger.  It should be at least three times 
#   the largest value observed in present-day species + 1. ", immediate.=TRUE)

  isNotNullGeomMean=!is.null(geomMean)
  isNotNullDirac=!is.null(dirac)
  if( (isNotNullGeomMean + isNotNullDirac + useRootStateMLE) !=1 )
    stop("Use exactly one of these three arguments: geomMean or dirac or useRootStateMLE")
  if (isNotNullDirac){
    isDiracInteger=(dirac==floor(dirac))
    if (!isDiracInteger){
      dirac<-floor(dirac)
      warning("the dirac value has to be an integer. It has been rounded down",immediate.=TRUE)
    }
    if (dirac<1) stop("the dirac value needs to be positive.")
  }
  if (isNotNullGeomMean){
    if (geomMean<1)
      stop("the mean of the prior geometric distribution has to be greater or equal to 1")
  }
  if (useRootStateMLE){
   if (conditioning=="oneInBothClades" | conditioning=="twoOrMore")
    stop(paste0("The 'rootStateMLE' is not implemented with the conditioning '",conditioning,"'"))
  }

  input = processInput(tr, equalBDrates, fixedRetentionRates, startingBDrates, startingQ)		
  
  nWGD=sum(input$phyloMat$type == "WGD")
  nWGT=sum(input$phyloMat$type == "WGT")
  nWG = nWGD + nWGT

  if ((nWG==0) & (!fixedRetentionRates)) {		
    fixedRetentionRates=TRUE
    warning("The tree has no WGD or WGT events, setting fixedRetentionRates to TRUE",
            immediate.=TRUE)
  } # important when the species tree has no WGD		
  para = input$para   # initial values for parameter
  lower = input$lower
  upper = input$upper
  if (!is.null(geomMean)){ geomProb=1/geomMean } else {geomProb=NULL} # processInput checked geomMean>=1.	
  conditioning=match.arg(conditioning)
  if (conditioning=="twoOrMore") {	
    geneNumberInEachFamily<-rowSums(geneCountData, dims = 1)
    familiesNotAppropriate=which(geneNumberInEachFamily==1)
    if (length(familiesNotAppropriate)>0){	
      print("The following families are responsible for the error:\n  ")
      print(familiesNotAppropriate)	
      stop("conditioning: we cannot use the type twoOrMore
 since at least one family has only 1 gene copy")
    }
  } else if (conditioning=="oneInBothClades") {
    .checkClades(input$phyloMat, geneCountData, input$nLeaf)
  }
  result <- optim(para, getLikGeneCount, input=input, geneCountData=geneCountData,
                  mMax=mMax, geomProb=geomProb, dirac=dirac, useRootStateMLE=useRootStateMLE,
                  conditioning=conditioning, fixedRetentionRates=fixedRetentionRates,
                  equalBDrates=equalBDrates, method="L-BFGS-B", lower=lower, upper=upper)	
  para = result$par
  loglik=-result$value

  if(equalBDrates){
    lambda = exp(para[1])
    lenlammu = 1 # length of vector for lambda and mu
    mu = "same as birth rate"
  } else { # not equal birth and deathrate
    lambda = exp(para[1])
    mu     = exp(para[2])
    lenlammu = 2
  }
  if( !fixedRetentionRates){ # update wgdTab with possibly new values in para
    input$wgdTab$retain2 = para[(lenlammu+1) : length(para)]
    input$wgdTab$retain1 = 1 - input$wgdTab$retain2
    if (nWGT){ # WGT: assuming independent loss of the 2 extra copies
      iT = input$wgdTab$type=="WGT"
      input$wgdTab$retain3[iT] =   input$wgdTab$retain2[iT]^2
      input$wgdTab$retain2[iT] = 2*input$wgdTab$retain2[iT]  * input$wgdTab$retain1[iT]
      input$wgdTab$retain1[iT] =   input$wgdTab$retain1[iT]^2
    }
  }
  return(list(birthrate=lambda, deathrate=mu, loglikelihood=loglik, WGDtable=input$wgdTab, 
              phyloMat=input$phyloMat, call=match.call(),convergence=result$convergence,mMax=mMax))
}
