processInput <- function(tr, equalBDrates=FALSE, fixedRetentionRates=TRUE,
                         startingBDrates=c(0.01, 0.02),startingQ=NULL)
{ # tr = phylo4d object from read.simmap()
  phyloMat = as.data.frame(tr@edge) 
  ## phyloMat is a phylo matrix with nNode rows and 6 columns: 
  ## Parent Child  Time Species RetenRate LossRate	 
  phyloMat = cbind(phyloMat, tr@edge.length)
  colnames(phyloMat) = c("Parent", "Child", "Time")	
  phyloMat = phyloMat[order(phyloMat$Child), ]
  ## order the phyloMat$Child column in an increasing order	
  phyloMat = cbind(phyloMat, Species=tr@label)
  r = which(phyloMat$Species == "Root") # root node
  phyloMat$Time[r] = -1                 # root has time -1

  y = as.character(tr@data[,1]) # y contains wgd vs. wgt info, but misses the root
  y = c(y[1:(r-1)],"rootPrior",y[r:length(y)])
  phyloMat$type = "BD"
  phyloMat$type[r] = "rootPrior"
  phyloMat$type[y=="wgd" | y=="WGD"] = "WGD"
  phyloMat$type[y=="wgt" | y=="WGT"] = "WGT"

#  retenRate = c(y[1:(r-1)], 0, y[r:length(y)])	
#  phyloMat = cbind(phyloMat, RetenRate=retenRate)
#  phyloMat$LossRate = 1 - phyloMat$RetenRate

  nLeaf = r-1
  nNode = nrow(phyloMat)
  nWGD = sum(phyloMat$type == "WGD")  # number of WGDuplications
  nWGT = sum(phyloMat$type == "WGT")  # number of WGTriplications
  nWG  = nWGD + nWGT
  iWG  = which(phyloMat$type =="WGD" | phyloMat$type =="WGT")
  wgdNode = phyloMat$Parent[iWG] # vector of "before WGD/T" nodes
  if (any(phyloMat$Time[iWG] != 0)){ 
    warning("One (or more) WGD or WGT edge had non-zero length,
  which will be set to 0.")
    phyloMat$Time[iWG] = 0
  }
  if (nWG==0) cat("found no WGD (or WGT)\n")
  else { cat("found: "); if (nWGD){
      cat(nWGD,"WGD(s)")
      if (nWGT) cat(" and",nWGT,"WGT(s).\n") else cat("\n")
    } else cat(nWGT,"WGT(s).\n")
  }

  ## "para": vector of starting values for later optimization,
  ## "lower" and "upper": bounds for "para"
  if (equalBDrates){
    para = log(startingBDrates[1])
    lower = c(-Inf)
    upper = c(Inf)
  } else {
    para = log(startingBDrates[1:2])
    lower = c(-Inf, -Inf)
    upper = c(Inf, Inf)
  }
  if (is.null(startingQ)) # setting all retention rates to 0.5
    startingQ = rep(0.5, nWG)
  else {                  # user-supplied retention rates
    if (length(startingQ)!= nWG)
      stop("Need ",nWG," retention rates in startingQ: number of WGD/Ts found in the tree.")
    if (nWG){
      if (sum(startingQ>1) > 0)
        stop("starting values for retention rates must be <= 1.")
      if (sum(startingQ<0) > 0)
        stop("starting values for retention rates must be >= 0.")
    }
  }
  if (!fixedRetentionRates){ # retention rate(s) to be included for optimization
    lower = c(lower, rep(0, nWG))
    upper = c(upper, rep(1, nWG))
    para = c(para,startingQ)
  }

  wgdTab = phyloMat[iWG, c("Parent","type"),drop=F] #, "RetenRate", "LossRate")]
  names(wgdTab)[1] = "wgdNode" # 1st column has node *before* the WGD/T
  wgdTab$retain1 = 1-startingQ # for WGDs: retain1 = loss rate
  wgdTab$retain2 =   startingQ #           retain2 = retention rate
  wgdTab$retain3 = rep(0,nWG)
  iT = wgdTab$type=="WGT"
  if (!is.null(iT)){
    wgdTab$retain1[iT] = (1-startingQ[iT])^2 # for WGTs: assumes independent loss of
    wgdTab$retain2[iT] = 2*startingQ[iT]*(1-startingQ[iT]) # the 2 extra copies
    wgdTab$retain3[iT] =   startingQ[iT]^2
  }
  
  edgeOrder = getEdgeOrder(phyloMat, nLeaf, wgdTab) # for post-order traversal
  return( list(phyloMat=phyloMat, nLeaf=nLeaf, nNode=nNode, wgdTab=wgdTab, edgeOrder=edgeOrder,
               para=para, lower=lower, upper=upper) )
}

getEdgeOrder <- function(phyloMat,nLeaf,wgdTab){
 # output: table with edges ordered for a post-order traversal,
 #         with info on which edges the birth-death process applies
 #         or a whole genome duplication/triplication event.
 # !! Warning !! assumes that:
 # 1. speciation nodes are given lower indices than
 #    singleton nodes when the tree is read in by phyext,
 # 2. *speciation* nodes are in post-order when going down from
 #    node maxInteNoSingNode to the root, whose index is nLeaf+1.
 # 3. the user correctly coded the nodes with WGD/T events,
 #    i.e. 2 singleton nodes for each WGD/T. 
  nNode = max(phyloMat$Child)	
  root = nLeaf +1
  maxInteNoSingNode = nNode - 2*(nrow(wgdTab))
  orderedNode = rep(NA,nNode)
  orderedEdge = rep(NA,nNode)     # index of parent edge in phyloMat
  ordered2sib = rep(NA,nNode)     # is node second sibling?
  orderedEdgeType = rep(NA,nNode) # WGD vs BD=birth-death vs. edge above root

  j=1 # next item in the ordered list
  for (parent in maxInteNoSingNode : root) { 
    parentIndex = which(phyloMat$Parent == parent)
    ## since parent is in inte-NoSing-Node, length(parentIndex) is 2, i.e. 2 children
    edge1 = parentIndex[1]; edge2 = parentIndex[2]
    child1 = phyloMat$Child[edge1]
    child2 = phyloMat$Child[edge2]
    branchNode1 = child1
    branchNode2 = child2
    # with what comes next, each branchNode will be a vector of nodes along
    # one "true" branch, including singleton nodes for WGD events, but excluding 'parent'. 
    # ex: if no singleton nodes on the branch, length(branchNode) = 1
    #     if 4 singletons (=2 WGD) on the branch, length(branchNode) = 5
    while(child1 > maxInteNoSingNode){
      nextedge = which(phyloMat$Parent==child1)
      child1 = phyloMat$Child[nextedge]
      branchNode1 = c(child1,branchNode1)
      edge1 = c(nextedge,edge1)
    }
    while(child2 > maxInteNoSingNode){
      nextedge = which(phyloMat$Parent==child2)
      child2 = phyloMat$Child[nextedge]
      branchNode2 = c(child2,branchNode2)
      edge2 = c(nextedge,edge2)
    }
    if (child1<=nLeaf){nc1=0} else {nc1=2} # again assuming binary tree
    if (child2<=nLeaf){nc2=0} else {nc2=2}
    l1 = length(branchNode1); l2=length(branchNode2)
    if (l1%%2==0 | l2%%2==0){
      stop("error in reading the species with the number of edges associated with WGDs.")}
    endj = j+l1+l2-1
    orderedNode[j:endj] = c(branchNode1,branchNode2)
    orderedEdge[j:endj] = c(edge1,edge2)
    ordered2sib[j:endj] = c(rep(FALSE,l1+l2-1),TRUE); 
    orderedEdgeType[j:endj] <- phyloMat$type[c(edge1,edge2)]
    ## c(rep(c("BD","WGD"),(l1-1)/2),"BD", rep(c("BD","WGD"),(l2-1)/2),"BD")
    j=endj+1
  }
  if (j!=nNode){ stop("error in likelihood calculation in post-order tree traversal")}
  orderedNode[j] = root # adding the root
  orderedEdge[j] = which(phyloMat$Child==root)
  ordered2sib[j] = FALSE # useless?
  orderedEdgeType[j] = "rootPrior"
  return(data.frame(child=orderedNode,edge=orderedEdge,type=orderedEdgeType,scdsib=ordered2sib))
}
