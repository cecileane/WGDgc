.BDprob <-
function (logLamlogMu, ti, mStart, nEnd) {
	# calculate probability from mStart genes in a parent to nEnd genes in a child in time=ti
        # birth death model
	
	if (mStart==0 && nEnd==0) {return(1)}
	
	if ( (length(logLamlogMu)==1) || (logLamlogMu[1] == logLamlogMu[2]) ){
		return( .calcProbSamePara(logLamlogMu[1], ti, mStart, nEnd) )
	}
	
	logLam = logLamlogMu[1]; logMu = logLamlogMu[2];
	
	lmdiff = exp(logLam) - exp(logMu)	# lam - mu
	mlratio = exp(logMu) / exp(logLam)	# mu/lam
	
	beta = (exp(lmdiff *ti) -1) / (exp(lmdiff *ti) -mlratio)
	alpha = beta * mlratio
	
	minMN = min(mStart,nEnd)
	value = choose(mStart+nEnd-1, mStart-1) *alpha^mStart *beta^nEnd 	# j=0 in the formula
	Pmn = value
	ab = (1-alpha-beta)/(alpha*beta)
	
	if (minMN >=1) {
		for (j in 1:minMN) {
			value = value *(mStart-j+1)/j *(nEnd-j+1)/(mStart+nEnd-j) *ab
			Pmn = Pmn + value
		}
	}
	return (Pmn)
}
.calcProbOneBranch <-
function(branchNode, logLamlogMu, nPos, nFamily, nLeaf, 
	geneCountData, wgdTab, phyloMat, MatDoomed, MDcol, Mat){
		
	 # branchNode is a vector of nodes on a branch in correct order as on the tree
	 #  ex:(6,15,14,13,12,7) denotes 6 is parent of 15, 15 is parent of 14 ... 12 is parent of 7
	 
	 # branchNode[1] = ancestor of this branch, it is internal-nonsingleton node
	 # branchNode[length(branchNode)] = last decendants of this branch, it is either 
	 #  a leaf or internal-nonsingleton node
	 
	 # if length(branchNode) >2, then this branch has some singletons in the middle
	 # MDcol is the column of MatDoomed, either 2 or 3
	 # this function returns an nPos x nFamily matrix of prob of data below ancestor node of this branch
 
	for(i in (length(branchNode) -1) :1){
		 # assume the prob of data below the last decendant is available at this point
		 # so calculate from the node next to the last up to the oldest ancestor
		
		parent = branchNode[i]
		child = branchNode[i+1]
		time = phyloMat$Time[phyloMat$Child==child]
		prob = matrix(0, nrow=nPos, ncol=nFamily)
		
		if (child <= nLeaf) { # this child is a leaf
			species = phyloMat$Species [phyloMat$Child == child]
			prob = .getPt(logLamlogMu=logLamlogMu, ti=time, nPos=nPos, nFamily=nFamily, isChildLeaf=T,
					nLeafCount=geneCountData[ ,species])

			MatDoomed[parent,MDcol]=prob[2,nFamily]

			if (length(branchNode)>2) {				
				#the parent is not an internal-nonSingleton node
				MatDoomed[parent,ifelse(MDcol==2,3,2)]=1
				MatDoomed[parent,1]=MatDoomed[parent,MDcol] 			
			}
			
			# getPt = function (logLamlogMu, ti, nPos, nFamily=NULL, isChildLeaf=F, nLeafCount=NULL, 
			#   isWgdNode=F, wgdLR=NULL) {
						
		} else { # child is an internal node
			if (time == 0) { # this parent is the node before wgd
				wgdLR = wgdTab$LossRate[wgdTab$wgdNode==parent]
				
				Pt = .getPt (logLamlogMu=logLamlogMu, ti=time, nPos=nPos, isWgdNode=T, wgdLR=wgdLR)
				MatDoomed[parent,ifelse(MDcol==2,3,2)]=1
				MatDoomed[parent,MDcol]= wgdLR*MatDoomed[child,1] +(1-wgdLR)*MatDoomed[child,1]^2
				MatDoomed[parent,1]=MatDoomed[parent,MDcol]
				
				
			} else { # parent is not node before WGD and child is not a leaf
				Pt = .getPt (logLamlogMu=logLamlogMu, ti=time, nPos=nPos)
				MatDoomed[parent,MDcol]=Pt[2,1]
				
				for (j in 2:nPos) {
				 	MatDoomed[parent,MDcol]=MatDoomed[parent,MDcol]+ Pt[2,j] * MatDoomed[child,1]^(j-1)
				}	
				
				if(i!=1){
					MatDoomed[parent,ifelse(MDcol==2,3,2)]=1
					MatDoomed[parent,1]=MatDoomed[parent,MDcol]
					}

			}
				
			for (k in 1:nPos) {
				prob[k, ] = Pt[k, ] %*% Mat[ ,child-nLeaf, ]
			}
		}
		
		Mat[ ,parent-nLeaf, ] = prob
	}
	return(list(prob=prob, MatDoomed=MatDoomed, Mat=Mat))
}
.calcProbSamePara <-
function (logLam, ti, mStart, nEnd) {
	# calculate probability from mStart genes in a parent to nEnd genes in a child in time=ti
	# here lamda = birth rate = death rate
	
	if (mStart==0 && nEnd==0) {return(1)}
	
	alpha = (exp(logLam) *ti) / (1+ exp(logLam)*ti)
	minMN = min(mStart, nEnd)
	value = choose(mStart+nEnd-1, mStart-1) *alpha^(mStart+nEnd)	# j=0 in the formula
	Pmn = value
	ab = (1-2*alpha) / alpha^2
	
	if (minMN >= 1) {
		for (j in 1:minMN){
			value = value *(mStart-j+1)/j *(nEnd-j+1)/(mStart+nEnd-j) *ab
			Pmn = Pmn + value
		}
	}
	return (Pmn)
}



.getPt <-
function (logLamlogMu, ti, nPos, nFamily=NULL, isChildLeaf=F, nLeafCount=NULL, isWgdNode=F, wgdLR=NULL) {
	# nLeafCount is 1 column in the geneCountData, it is a vector of genecounts of 1 species in all families
	# nPos = max posible number of genes at each internal node, for ex: nPos=101 (0 to 100 genes) 
	#   recommend choosing nPos=2*max(geneCountData)
		
	# this function return Pt which is a probability matrix
	
	# if isChildLeaf=T, Pt is a nPos x nFamily matrix in which each entry Pt[i,j]=prob from i-1 genes 
	#   to nLeafCount[j] genes in family j 
	# if isChildLeaf=F, Pt is nPos x nPos matrix, each entry Pt[i,j] = prob from i-1 genes to j-1 genes
	#   (based on birth-death process)
	
	# if isWgdNode=T, Pt is nPos x nPos matrix, each entry Pt[i,j] = prob from i-1 genes to j-1 genes
	#   (based on binomial distribution)
	# isWgdNode = whether or not this is a node before WGD event
	# wgdLR = wgd loss rate
		
	
	if (isWgdNode) {
		Pt = matrix (0, nrow=nPos, ncol=nPos) 
		for (i in 0:(nPos-1)) {
			N = ifelse( 2*i > (nPos-1), nPos-1, 2*i)
			
			for (j in i:N) {
				Pt[i+1,j+1] = choose(i, j-i) *(1-wgdLR)^(j-i) *wgdLR^(2*i -j)
				 # binomial dist of getting j genes (after wgd) from i genes (before wgd),
				 # j is between i:2i or i:(nPos-1)
			}
		}
		return (Pt)
	}	
	
	if (isChildLeaf) {
		Pt = matrix(0, nrow=nPos, ncol=nFamily)
		for (i in 1:nPos) {
			for (j in 1:nFamily) {
				Pt[i,j] = .calcProb (logLamlogMu, ti, i-1, nLeafCount[j])
			}
		}
	} else {
		Pt = matrix (0, nPos, nPos)
		for (i in 1:nPos) {
			for (j in 1:nPos) {
				Pt[i,j] = .calcProb (logLamlogMu, ti, i-1, j-1)
			}
		}
	}
	# calcProb = function (logLamlogMu, ti, mStart, nEnd) {}
	return (Pt)
}


.checkClades <-
function(phyloMat, geneCountData, nLeaf){
	#check if all families have at least 1 gene copy in each clade

	r = which(phyloMat$Species == "Root")
	descRoot=which(phyloMat$Parent==r)
	tableClade<-array (0, nLeaf)

	for (Leafnode in 1 : nLeaf) {

		node<-Leafnode
		
		if (node==descRoot[1]){
        	tableClade[Leafnode] =1 
		} else if (node==descRoot[2]) {
       	 tableClade[Leafnode] = -1
		}

		while (node!=descRoot[1] & node!=descRoot[2]) { 
	
			indnode<-which(phyloMat$Child==node)
			node<-phyloMat$Parent[indnode]
			if (node==descRoot[1]){
        			tableClade[Leafnode] =1 
			} else if (node==descRoot[2]) {
				tableClade[Leafnode] = -1
			}
		}	
	}


	codeClade=-1
	#we begin with the first clade

	for (clade in 1:2){
	#loop on the different clades

		indClade=which(tableClade==codeClade)
		SpecName<-phyloMat$Species[indClade]
		ind<-which(colnames(geneCountData)%in%SpecName)
		geneCountDataClade <- geneCountData[, c(ind)]

		if (length(ind)>1) {
			geneNumberInEachFamilyClade<-rowSums(geneCountDataClade, dims = 1)
		}else {
			geneNumberInEachFamilyClade<-geneCountDataClade
		}

		indEmptyFamily=which(geneNumberInEachFamilyClade==0)

		if (length(indEmptyFamily)>0){
			print("The following families are responsible for the ERROR")
			print(indEmptyFamily)
			stop("ERROR about the conditioning, we can not use the type oneInBothClades 
		     	since at least one family has 0 gene copies in one of the 2 clades")

		}else{
			codeClade=1
			#we have to check the second clade
		}

	}

}
