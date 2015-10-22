logLik_CsurosMiklos <-
function(logLamlogMu, nLeaf, nFamily, phyloMat, geneCountData, mMax, wgdTab, edgeOrder) {
 # by Cecile Ane, October 2013. Based on Csuros & Miklos (2009).
 # Replaces getMatAndMatDoomed(). Used by getLikGeneCount().
 # phyloMat, wgdTab: see processInput.r. edgeOrder: see getEdgeOrder().

 nNode = max(phyloMat$Child)	
 root = nLeaf +1

 sameBDrates = length(logLamlogMu)==1 || isTRUE(all.equal(logLamlogMu[1],logLamlogMu[2]))
 logLam = logLamlogMu[1]
 #cat("  lambda=",exp(logLam),"\n")
 if (!sameBDrates){
  logMu = logLamlogMu[2]
  #cat("  mu    =",exp(logMu),"\n")
  lmdiff = exp(logLam) - exp(logMu)	# lam - mu
  mlratio = exp(logMu) / exp(logLam)	# mu/lam
 }

 doomedNode = array (NA, nNode)
 # probability that a lineage starting at the node is doomed below that node: D values
 doomedEdge = array (NA, nNode)
 # probability that a lineage starting at the beginning of the parent edge
 # of the node is doommed below that node. G_xi(0) values in Csuros & Miklos (2009).

 transitionProb = array(0, c(mMax+1,mMax+1, nNode))
 # entry (i+1,j+1,n) is the probability that there are j genes at node n
 # and that the i original genes all survive below node n, given that 
 # there are i genes at the parent of node n. w* values in Csuros & Miklos (2009)
 transitionProb[1,1, ]= 1 # w*(0|0)=1 for all nodes. Others initialized to 0


 loglik = array (-Inf, c(mMax+1, nNode-nLeaf, nFamily))
 # 3D matrix of size M x (nNode-nLeaf) x nFamily
 # entry (i,j,k): probability of data below node j+nLeaf in family k
 # given that there are i surviving genes at node j+nLeaf
 # !!Warning: assumes 2 children at each speciation node in the species tree!!

 # Auxiliary value B_n: probability of observing data below node n and
 # first s genes at the parent of n all survive below n, given that
 # there are s+t genes at the parent of n. Not needed on WGD edges.
 # case t=0 only needed for first sibling edges. 
 # 2 data structures to save space:
 # logB0 = log of B_node for t=0 for all edges, and
 # logB for t>0 only and for second siblings only.
 sib2toNode = edgeOrder$child[ edgeOrder$scdsib]
 nsib2 = length(sib2toNode)
 node2sib = rep(NA,nNode)
 node2sib[sib2toNode] = 1:nsib2
 logB  = array(NA, c(mMax+1, nsib2, nFamily, mMax))
 # log(B): s=i-1, node=j, family=k, t=ell, i.e. starts at t=1
 logB0  = array(NA, c(mMax+1, nNode, nFamily))

 for (ie in 1:dim(edgeOrder)[1]){
  
  node = edgeOrder$child[ie]
  edge = edgeOrder$edge[ie]
  childEdges = which(phyloMat$Parent == node)
  childNodes = phyloMat$Child[childEdges]

  # calculate doomedNode: Dx=prod Gx of children
  if (length(childEdges)){ # node x is not a leaf
   doomedNode[node] = prod( doomedEdge[childEdges] )
  } else {                # node x is a leaf
   doomedNode[node]=0 # Change this to account for data error model
  }

  # calculate lokLik: from B values of children edges
  if (length(childEdges)==1){   # node not a leaf, but single child edge
   loglik[,node-nLeaf,] = -(0:mMax)*log(1-doomedEdge[childEdges]) + logB0[,childNodes,]
   #                        0:mMax correctly recycled to fill in array
  }
  if (length(childEdges)>1) { # there should be 2 children
   # first part: sum n!/(s!(n-s)!) D_x1^s B_x1;0,n-s  B_x2;n-s,s
   logDedge1 = log(doomedEdge[childEdges[1]])
   # lik: mMax+1 x nFamily matrix, node specific.
   # initializing with term from t=0=n-s: D_x1^s B_x1;0,0  B_x2;0,s
   lik = exp((0:mMax)*logDedge1 +rep(logB0[0+1,childNodes[1],],each=mMax+1) +logB0[,childNodes[2],])
   # adding contribution of all t=n-s>0 values, for s fixed.
   for (s in 0:(mMax-1)){ # s=mMax contributes to n=mMax only, with t=0: already covered
    tMax = mMax-s
    lik[s +1:tMax +1,] = lik[s +1:tMax +1,] +
      exp( lchoose(1:tMax+s,s) + s*logDedge1 + logB0[1:tMax +1,childNodes[1],] + 
           t(logB[s+1,node2sib[childNodes[2]],,1:tMax]))
   }
   # second part:  * (1-Dx)^(-n)
   loglik[,node-nLeaf,] = log(lik) - (0:mMax)*log(1-doomedNode[node])
  }

  # calculate doomedEdge and transition probs on parent edge.
  if (edgeOrder$type[ie] == "BD"){
   len = phyloMat$Time[edge] # branch length
   if (sameBDrates){
    gam = exp(logLam)*len/(1+exp(logLam)*len)
    psi   = gam
   } else {
    psi = (exp(lmdiff *len) -1) / (exp(lmdiff *len) -mlratio)
    gam = psi * mlratio
   }
   doomedEdge[edge] = gam + (1-gam)*(1-psi)*doomedNode[node]/(1-psi*doomedNode[node])
   psiprime = 1 - (1-psi)/(1-psi*doomedNode[node])
   G1 = (1-doomedEdge[edge])*(1-psiprime)
   for (j in 1:mMax){ # w*(j|i) = G1 w*(j-1|i-1) + psi' w*(j-1|i), for j>=i.
    transitionProb[1:j +1, j+1, node] = G1* transitionProb[1:j,  j, node] +
                                  psiprime* transitionProb[1:j+1,j, node]
   }
  }
  if (edgeOrder$type[ie] == "WGD"){
   ievent = wgdTab$wgdNode == phyloMat$Parent[edge]
   p1g = wgdTab$retain1[ievent] # prob to retain 1 gene
   p2g = wgdTab$retain2[ievent]
   doomedEdge[edge] = p1g * doomedNode[node] + p2g * doomedNode[node]^2
   G1 = (p1g + 2*p2g*doomedNode[node]) * (1-doomedNode[node])
   G2 = p2g * (1-doomedNode[node])^2
   transitionProb[2,2, node]=G1 # w*(1|1)
   transitionProb[2,3, node]=G2 # w*(2|1)
   # transitionProb[1+ 0:mMax * (mMax+2) + (mMax+1)*(mMax+1)*(node-1)] <- G1^(0:mMax)
   for (i in 2:mMax){ # w*(j|i) = G1 w*(j-1|i-1) + G2 w*(j-2|i-1), for j in [i,2i].
    jmax = min(2*i, mMax)
    transitionProb[i+1, i:jmax+1, node] = G1* transitionProb[i, i:jmax,   node] +
                                          G2* transitionProb[i, i:jmax-1, node]
   }
  }
  if (edgeOrder$type[ie] == "WGT"){
   ievent = wgdTab$wgdNode == phyloMat$Parent[edge]
   p1g = wgdTab$retain1[ievent]
   p2g = wgdTab$retain2[ievent]
   p3g = wgdTab$retain3[ievent]
   doomedEdge[edge] = p1g * doomedNode[node] + p2g * doomedNode[node]^2 + p3g * doomedNode[node]^3
   G1 = (p1g + 2*p2g*doomedNode[node] + 3*p3g*doomedNode[node]^2) * (1-doomedNode[node])
   G2 = (p2g + 3*p3g*doomedNode[node]) * (1-doomedNode[node])^2
   G3 = p3g * (1-doomedNode[node])^3
   transitionProb[2,2, node]=G1 # w*(1|1)
   transitionProb[2,3, node]=G2 # w*(2|1)
   transitionProb[2,4, node]=G3 # w*(3|1)
   # transitionProb[1+ 0:mMax * (mMax+2) + (mMax+1)*(mMax+1)*(node-1)] <- G1^(0:mMax)
   for (i in 2:mMax){ # w*(j|i) = G1 w*(j-1|i-1) + G2 w*(j-2|i-1) + G3 w*(j-3|i-1), for j in [i,3i].
    jmax = min(3*i, mMax)
    if (i>2){
     transitionProb[i+1, i:jmax+1, node] = G1* transitionProb[i, i:jmax,   node] +
                                           G2* transitionProb[i, i:jmax-1, node] +
                                           G3* transitionProb[i, i:jmax-2, node]
    } else { # i=2. problem for j=i: j-3 is outside the bounds of table
     transitionProb[2+1, 2+1, node] = G1* transitionProb[2, 2,   node] +
                                      G2* transitionProb[2, 2-1, node]
     if (jmax >=3)
       transitionProb[2+1,3:jmax+1,node] = G1* transitionProb[2, 3:jmax,   node] +
                                           G2* transitionProb[2, 3:jmax-1, node] +
                                           G3* transitionProb[2, 3:jmax-2, node]
    }
   }
  }

  # calculate Bs along parent edge, unless root edge
  if (edgeOrder$type[ie] != "rootPrior"){ # B values for t=0
   if (length(childEdges)==0){ # leaf edge: B0 = transitionProb(s,data at the tip)
    leafName = as.character(phyloMat$Species[edge])
    logB0[,node,] = log( transitionProb[,geneCountData[,leafName]+1,node] )
    fam.na = which(is.na(geneCountData[,leafName])) # families with missing data at the leaf
    logB0[,node,fam.na] = (0:mMax)*log(1-doomedEdge[edge]) # 0:mMax recycled okay
   } else {
    logB0[,node,] = log( transitionProb[,,node] %*% exp(loglik[,node-nLeaf,]) )
   }
   if (edgeOrder$scdsib[ie]){   # other B values, for t>0
    lde = log(doomedEdge[edge]) # log of Gx(0)
    # initialize with case s=mMax: B_t,mMax = Gx(0)^t * B_0,mMax
    logB[mMax+1,node2sib[node],,1:mMax] =
      matrix( (1:mMax)*lde, nFamily,mMax,byrow=T)  +
      logB0[mMax+1,node,] # logB0 of families recycled okay across m values
    for (t in 1:mMax){ # B_t,s = B_t-1,s+1 + Gx(0)*B_t-1,s
     #logB[1:mMax,node,,t+1] = log(exp(logB[(1:mMax)+1,node,,t]) + exp(logB[1:mMax,node,,t]+log(doomedEdge[edge])))
     # to limit rounding errors, calculate log( exp(x) + exp(y) ) as: max + log(1 + exp(min-max))
     # but issues when x=-Inf and y=-Inf, because min-max = -Inf - -Inf = NaN
     if (t>1){
      x=logB[(1:mMax)+1,node2sib[node],,t-1]
      y=logB[ 1:mMax,   node2sib[node],,t-1] + lde
     } else {
      x=logB0[(1:mMax)+1,node,]
      y=logB0[ 1:mMax,   node,] + lde
     }
     tmp2 = pmax(x, y)
     tmp3 = which(tmp2 != -Inf)
     logB[1:mMax,node2sib[node],,t][tmp3] = tmp2[tmp3] + log(1+exp(pmin(x[tmp3], y[tmp3]) - tmp2[tmp3]))
     logB[1:mMax,node2sib[node],,t][-tmp3]= -Inf
    }
   }
  }

 }
 #return(list(loglik=loglik,doomedNode=doomedNode,doomedEdge=doomedEdge,w=transitionProb,logB0=logB0,logB=logB))
 return(list(loglikRoot=loglik[,root-nLeaf,],doomedRoot=doomedNode[root],
             doomedRootLeft =doomedEdge[childEdges[1]],doomedRootRight=doomedEdge[childEdges[2]]))
}

# .logsumexp <- function(x,y){ # gets log( exp(x) + exp(y) ) carefully
#  tmp = pmax(x,y); return(tmp + log(1+exp(pmin(x,y)-tmp))) }
