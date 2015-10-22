.logLikGeneCount <-
function (para, input, geneCountData, nPos=NULL, geomMean=NULL, dirac=NULL, useRootStateMLE=F, conditioning=c("oneOrMore", "twoOrMore", "oneInBothClades", "none"), equalBDrates=F, fixedRetentionRates=T) {
	# para is a vector of parameters
	#   para=c(logLam) if equalBDrates=T and fixedRetentionRates=T
	#   para=c(logLam, logMu) if equalBDrates=F and fixedRetentionRates=T
	#   para=c(logLam, q1, q2,...) if equalBDrates=T and fixedRetentionRates=F
	#   para=c(logLam,logMu, q1, q2,...) if equalBDrates=F and fixedRetentionRates=F
	
	# geneCountData is a dataframe which each row corresponds to a gene family and each column corresponds to a species
	# geomMean is the mean of the prior geometric dist, dirac is the mean of prior dirac dist
	# useRootStateMLE = whether or not choosing the best number of gene to start with
	
	#conditioning: 
	#conditioning="oneOrMore" if all families with one or more gene copies are included in the data, 
	#conditioning="twoOrMore" to condition on families having two of more genes, 
	#conditioning="oneInBothClades" if the data set was filtered to include only families with at least one gene copy 
	#in each of the two main clades stemming from the root.
 	#conditioning="none" uses unconditional likelihoods.
	
	# input: read processInput.R for explanation
	# this function return -log(likelihood)
	
	if(equalBDrates){
		logLamlogMu = para[1]
	} else {
		logLamlogMu = para[1:2]
	}
	
	wgdTab = input$wgdTab
	 # this table has all the nodes before WGD with their corresponding LossRate and RetenRate
	 
	if( !fixedRetentionRates){
		wgdTab$RetenRate = para[(length(logLamlogMu)+1) : length(para)]
		wgdTab$LossRate = 1 - wgdTab$RetenRate
	}

	nLeaf = input$nLeaf
	phyloMat = input$phyloMat
	
	if (conditioning=="oneOrMore") {
		geneCountData = rbind(geneCountData, rep(0, nLeaf))
		# add 1 row of zero's to geneCountData to calculate prob of observe nothing
		
	} else if (conditioning=="twoOrMore") {
		addiFamily = data.frame( rbind( diag(1, nrow=nLeaf, ncol=nLeaf), rep(0, nLeaf)))
		names(addiFamily) = names(geneCountData)
		# add nLeaf+1 additional families to the geneCountData
			
		geneCountData = rbind(geneCountData, addiFamily)
			
	} else if (conditioning=="oneInBothClades") {
		geneCountData = rbind(geneCountData, rep(0, nLeaf))
		# add 1 row of zero's to geneCountData to calculate prob of observing nothing
			
	} 
		
	nFamily = nrow(geneCountData)
		
	myResults = .getMatAndMatDoomed(logLamlogMu, nLeaf, nFamily, phyloMat, geneCountData, nPos, wgdTab)
	
	Mat=myResults$first
	MatDoomed=myResults$second

	if (!is.null(dirac)) {
		lkli = Mat[dirac+1, 1,]	
	} else if (useRootStateMLE) {
		lkli = apply( Mat[ 2:nPos,1, ], 2, max)
		 #print (apply( Mat[ 2:nPos,1, ], 2, function(x) { which( x==max(x) ) })) 	
	} else if (!is.null(geomMean)){
		rootPriorPmf = dgeom(0:(nPos-2), ifelse(geomMean >=1, 1/geomMean, .99))	
		 # prior prob mass at the root
		
		lkli = rootPriorPmf %*% Mat[ 2:nPos, 1, ]
		 # column 1 is the root node of the tree
		 # starting from row 2 since at the root, number of starting genes should be nonzero
	}
	 # lkli is a vector whose length equals the number of gene families, 
	 # each entry lkli[i]=likelihood of the the tree in family i
	 
	if (conditioning=="none") {
		logLkli = sum(log( lkli ))
	
	} else if (conditioning=="oneOrMore") {
		logLkli = sum(log( lkli[1:(nFamily-1)] )) - (nFamily-1)*log(1-lkli[nFamily])
		 # lkli[nFamily] is prob of observing nothing
	
	} else if (conditioning=="twoOrMore"){
		logLkli = sum(log( lkli[1:(nFamily-nLeaf-1)] )) - (nFamily-nLeaf-1)*log(1 - sum(lkli[ (nFamily-nLeaf):nFamily ]))
		 # lkli[ nFamily-nLeaf : nFamily] are prob of added families that have only 0 or 1 gene in total
	
	} else if (conditioning=="oneInBothClades"){		
		#denom = 1-(MatDoomed[nLeaf+1,2])/(1-MatDoomed[nLeaf+1,2])-(MatDoomed[nLeaf+1,3])/(1-MatDoomed[nLeaf+1,3])	+(MatDoomed[nLeaf+1,1])/(1-MatDoomed[nLeaf+1,1])
	 	
	 	denom = 1-MatDoomed[nLeaf+1,2]-MatDoomed[nLeaf+1,3]+MatDoomed[nLeaf+1,1]
		logLkli = sum(log( lkli[1:(nFamily-1)] ))  - (nFamily-1)*log(denom)
	}
		
	return (-logLkli) 
}

