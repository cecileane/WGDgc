getLikGeneCount <- function (para, input, geneCountData, mMax=NULL,
                             geomProb=NULL, dirac=NULL, useRootStateMLE=FALSE,
                             conditioning=c("oneOrMore", "twoOrMore", "oneInBothClades", "none"),
                             equalBDrates=FALSE, fixedRetentionRates=TRUE) {
  # para is a vector of parameters
  #   para=c(logLam) if equalBDrates=T and fixedRetentionRates=T
  #   para=c(logLam, logMu) if equalBDrates=F and fixedRetentionRates=T
  #   para=c(logLam, q1, q2,...) if equalBDrates=T and fixedRetentionRates=F
  #   para=c(logLam,logMu, q1, q2,...) if equalBDrates=F and fixedRetentionRates=F
  
  # geneCountData is a dataframe which each row corresponds to a gene family and each column corresponds to a species
  # 1/geomProb is the mean of the prior geometric dist, dirac is the mean of prior dirac dist
  # useRootStateMLE = whether or not choosing the best number of gene to start with
  
  #conditioning: 
  #conditioning="oneOrMore" if all families with one or more gene copies are included in the data, 
  #conditioning="twoOrMore" to condition on families having two of more genes, 
  #conditioning="oneInBothClades" if the data set was filtered to include only families with at least one gene copy 
  #in each of the two main clades stemming from the root.
  #conditioning="none" uses unconditional likelihoods.
  
  # input: read processInput.R for explanation
  # this function return -log(likelihood)

  if(equalBDrates)
    logLamlogMu = para[1]
  else 
    logLamlogMu = para[1:2]
 
  if( !fixedRetentionRates){ # update wgdTab with possibly new values in para
    input$wgdTab$retain2 = para[(length(logLamlogMu)+1) : length(para)]
    input$wgdTab$retain1 = 1 - input$wgdTab$retain2
    iT = input$wgdTab$type=="WGT"
    if (!is.null(iT)){ # WGT: assuming independent loss of the 2 extra copies
      input$wgdTab$retain3[iT] =   input$wgdTab$retain2[iT]^2
      input$wgdTab$retain2[iT] = 2*input$wgdTab$retain2[iT]  * input$wgdTab$retain1[iT]
      input$wgdTab$retain1[iT] =   input$wgdTab$retain1[iT]^2
    }
  }

  nLeaf = input$nLeaf
	
  if (conditioning=="twoOrMore") {
    addiFamily = data.frame( diag(1, nrow=nLeaf, ncol=nLeaf) )
    names(addiFamily) = names(geneCountData)
    ## add nLeaf additional families to the geneCountData	
    geneCountData = rbind(geneCountData, addiFamily)		
  }		
  nFamily = nrow(geneCountData)
        
  ##myResults = getMatAndMatDoomed(logLamlogMu, nLeaf, nFamily, phyloMat, geneCountData, nPos, wgdTab)
  myResults <- logLik_CsurosMiklos(logLamlogMu, nLeaf, nFamily, input$phyloMat,
                                   geneCountData, mMax, input$wgdTab, input$edgeOrder)
	
  ## Integrate over the prior at the root, to get log-likelihood of each gene family.
  if (!is.null(dirac)) {
    # lkli = Mat[dirac+1, 1,]
    imax = min(dirac,mMax) # max #of surviving lineages with >0 prior prob and >0 likelihood
    rootprior = dbinom(0:imax,size=dirac,prob=1-myResults$doomedRoot)
    loglik = log(rootprior %*% exp(myResults$loglikRoot[0:imax+1,]))
    doom = myResults$doomedRoot^dirac
    doomL = myResults$doomedRootLeft^dirac
    doomR = myResults$doomedRootRight^dirac
  } else if (useRootStateMLE) {
    ## lkli = apply( Mat[ 2:nPos+1,1, ], 2, max)
    loglik = apply(myResults$loglikRoot,2,max)
    doom <- 0 # assuming that there aren't any observed doomed families: with counts [0,...,0]
    ## fixit: for conditioning="oneInBothClades" or "twoOrMore"
    #         estimate the conditional likelihood here, and take the max *after* conditioning.
    # calculate doomL and doomR, for each estimated # of surviving lineages at the root.
    # doomL = prob{family doomed in left clade | n surviving genes at the root}
    #       = prob{gene counts = [0,...,0] on the left | n surviving genes at the root}
  } else if (!is.null(geomProb)){ # processInput already checked that geomProb>0 and <=1
    ## geomProb is for the prior number lineages at the root above.
    ## For the prior number of surviving lineages, modified shifted geometric distribution:
    prbprime = geomProb/(1-(1-geomProb)*myResults$doomedRoot) # this is 1-psiprime
    doom = prbprime * myResults$doomedRoot                    # this is gammaprime
    rootprior = c(doom , (1-doom)*dgeom(0:(mMax-1), prbprime))
    loglik = log(rootprior %*% exp(myResults$loglikRoot))
    doomL = geomProb * myResults$doomedRootLeft /(1-(1-geomProb)*myResults$doomedRootLeft)
    doomR = geomProb * myResults$doomedRootRight/(1-(1-geomProb)*myResults$doomedRootRight)
  }

  ## combine families and condition on data filtering	
  if (conditioning=="none") {
    loglik.all = sum(loglik)
  } else if (conditioning=="oneOrMore") {
    if (isTRUE(all.equal(doom,1)))
      stop("Doom probability is 1, cannot condition on non-extinct family.
  BD rates:",exp(logLamlogMu))  
    loglik.all = sum(loglik) - nFamily * log(1-doom)
  } else if (conditioning=="oneInBothClades"){
    ## Wrong in previous version: used the probability that a *single* lineage
    ## does not go extinct in both the left and right clade.
    ## need to account for prior distribution of # lineages at the root.
    if (isTRUE(all.equal(c(doomL,doomR),c(1,1))))
      stop("Doom probabilities are 1, cannot condition on non-extinct family in both clades.
  BD rates:",cat(exp(logLamlogMu)))
    loglik.all = sum(loglik) - nFamily * log(1-doomL-doomR+doom)
    # fixit: if (useRootStateMLE), doomL and doomR vary across families:
    #        depend on estimated # lineages surviving at root. See above.
  } else if (conditioning=="twoOrMore"){
    ## nFamily-nLeaf = real number of families. nLeaf families were added with 1 gene each.
    probObserved = 1 - doom - sum(exp(loglik[(nFamily-nLeaf+1):nFamily]))
    loglik.all = sum(loglik[1:(nFamily-nLeaf)]) - (nFamily-nLeaf)*log(probObserved)
  }	
  return (-loglik.all) 
}

