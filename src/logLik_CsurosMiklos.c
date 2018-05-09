#include "WGDgc.h"

/* fixit: change the function to access it in R with .Call instead of .C 
 see www.sfu.ca/~sblay/R-C-interface.ppt    
 SEXP getInt(SEXP myint, SEXP myintVar) {
  int Imyint, n; // declare an integer variable
  int *Pmyint;   // pointer to an integer vector
  PROTECT(myint = AS_INTEGER(myint));
  n = INTEGER_VALUE(myintVar);
  Pmychar[1] = R_alloc(strlen(CHAR(STRING_ELT(mychar, 1))),  sizeof(char)); 
  strcpy(Pmychar[0], CHAR(STRING_ELT(mychar, 0))); 
  return(R_NilValue);
  SEXP myint;
  int *p_myint;
   int len = 5;
 // Allocating storage space:
   PROTECT(myint = NEW_INTEGER(len));
   p_myint = INTEGER_POINTER(myint);
   p_myint[0] = 7;
   UNPROTECT(1);
   return myint;
 to work with real numbers, replace 
 int with double and INTEGER with NUMERIC
 compare to NA_INTEGER or NA_REAL or use macro ISNA
  book p.191: books.google.com/books?isbn=1420063685
*/

void logLik_CsurosMiklos (double *logLamlogMu, int *nLeaf, int * nFamily,
			  int *geneCountData, int *mMax,  int * equalBD,
			  int *nNode, double *edgelength, int *des, int *anc, int *scdsib, int *nsib2,
			  int *ievent, double *p1g, double *p2g, double *p3g) {
/* by Cecile Ane, October 2013 for the R version, June 2014 for C version (but this is not finished!!).
   Based on algorithm by Csuros & Miklos (2009).
   used by getLikGeneCount().
   ievent: vector of size nNode (one for each edge). nNode includes the root. 
           0 if "BD" type, -1 if root edge, R-index of WGD/WGD event otherwise.
   p1g, p2g, p3g: vectors of size # WGD or WGT events.
   nsib2: number of nodes who are second siblings.
   scdsib: 0 if the node is *not* the second sibling of its parent. R-index of first sibling otherwise
   Assumptions:
     speciation nodes have no more than 2 children edges.
     leaf IDs (R-style) are 1:nLeaf. Root ID = nLeaf+1.
     column i in geneCountData if for tip i (R-index)
     missing count data coded as -1.
*/
// fixit: extract edgelength, des, anc (ancestor), type, scdsib from phyloMat, wgdTab, edgeorder
// fixit: reorder genecountdata prior to this, so that data for des=i is in row i.

  int root = *nLeaf;

  // sameBDrates = length(logLamlogMu)==1 || isTRUE(all.equal(logLamlogMu[1],logLamlogMu[2]))
  // fixit: test equality in R, prior to this.
  double logLam = logLamlogMu[0];
  double logMu = logLam;
  double lmdiff, mlratio;
  if (! *equalBD){
    logMu = logLamlogMu[1];
    lmdiff = exp(logLam) - exp(logMu);	// lam - mu
    mlratio = exp(logMu) / exp(logLam);	// mu/lam
  }

  int nN = *nNode;
  double* doomedNode = (double*)calloc(nN, sizeof(double));
  // probability that a lineage starting at the node is doomed below that node: D values
  double* doomedEdge = (double*)calloc(nN, sizeof(double));
  // probability that a lineage starting at the beginning of the parent edge
  // of the node is doommed below that node. G_xi(0) values in Csuros & Miklos (2009).

  int mM = *mMax+1;
  int mM2 = mM * mM;
  int mM2nN = mM2 * nN;
  double* transitionProb = (double*)calloc( mM2nN, sizeof(double));
  /* 3-dimensional array.use C fr_finite
     entry (i,j,n+1) is the probability that there are j genes at node n
     and that the i original genes all survive below node n, 
     given i genes at the parent of n. w* values in Csuros & Miklos (2009) */
  for (int n=0; n<mM2nN; n+=mM2)
    transitionProb[k] = 1.0; // w*(0|0)=1 for all nodes.
  // Others initialized to 0 with calloc
  printf("transitionprob:\n")
  for (int n=0; n<mM2nN; n+=mM2){
    for (int j=0; j<mM2; j+=mM){
      for (int i=0; i<mM; i++)
        printf("%10d\n",transitionProb[i+j+n]);
        // transitionProb[i+j+k] = 0.0;
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n")
  double mInf = -1.0/0.0;
  int niN = nN - *nLeaf;
  int mMniN = mM * niN;
  int mMniNnF = mMniN * *nFamily;
  double* loglik = (double*)malloc( mMniNnF, sizeof(double));
  /* 3D matrix of size (mMax+1) x (nNode-nLeaf) x nFamily
     entry (i,j,k): probability of data below node j+nLeaf in family k
     given that there are i surviving genes at node j+nLeaf
     !!Warning: assumes 2 children at each speciation node in the species tree */
   // fixit: initialize this to -Inf
  for (int f=0; f<mMniNnF; f+=mMniN;){
    for (int n=0; n<mMniN; n+=mM;){
      for (int i=0; i<mM; i++){
        loglik[i+n+f] = mInf;
        printf("%10d\n", loglik[i+n+f]);
      }
      printf("\n");
    }
    printf("\n");
  }

  /* Auxiliary value B_n: probability of observing data below node n and
     first s genes at the parent of n all survive below n, given that
     there are s+t genes at the parent of n. Not needed on WGD edges.
     case t=0 only needed for first sibling edges. 
     2 data structures to save space:
     logB0 = log of B_node for t=0 for all edges, and
     logB for t>0 only and for second siblings only. */
  int mMnNnF = mMnN * nF;
  double* logB0 = (double*)calloc( mMnNnF, sizeof(double)); // initialized to 0

  int nS = *nsib2;
  int nF = *nFamily;
  int mMnS = mM * nS;
  int mMnSnF = mMnS * nF;
  int mMnSnFmMax = mMnSnF * *mMax; // !not! mMax+1 here.
  double* logB = (double*)calloc( mMnSnFmMax, sizeof(double));
  // log(B): s=i-1, node=j, family=k, t=ell, i.e. starts at t=1
  // sib2toNode = edgeOrder$child[ edgeOrder$scdsib]
  // node2sib = rep(NA,nNode)
  // node2sib[sib2toNode] = 1:nsib2

  for (int iedge=0; iedge<nN; iedge++){
    /* assume post-order traversal */
    int di = des[iedge]-1; // C-style index, of both the child node and the edge
    int ai = anc[iedge]-1; // for ancestor, or parent edge. fixit: use -1 for ancestor of root.
    // fixit: update the parent edge instead of searching for info from the child edges
    //childEdges = which(phyloMat$Parent == node)
    //childNodes = phyloMat$Child[childEdges]
    // add whatever to make sure loglik[,di-nLeaf,] and doomedNode[di] are correct at this point:
    // contributions from both children have been taken in and updated. 

    /* calculate doomedNode: Dx=prod Gx of children */
  if (length(childEdges)){ # node x is not a leaf
   doomedNode[node] = prod( doomedEdge[childEdges] )
  } else {                # node x is a leaf
   doomedNode[node]=0 # Change this to account for data error model
  }

  /* calculate lokLik: from B values of children edges */
  if (length(childEdges)==1){   # node not a leaf, but single child edge
   loglik[,node-nLeaf,] = -(0:mMax)*log(1-doomedEdge[childEdges]) + logB0[,childNodes,]
     //                        0:mMax correctly recycled to fill in array
  }
  if (length(childEdges)>1) { // there should be 2 children
    // first part: sum n!/(s!(n-s)!) D_x1^s B_x1;0,n-s  B_x2;n-s,s
   logDedge1 = log(doomedEdge[childEdges[1]])
     // lik: mMax+1 x nFamily matrix, node specific.
     // initializing with term from t=0=n-s: D_x1^s B_x1;0,0  B_x2;0,s
   lik = exp((0:mMax)*logDedge1 +rep(logB0[0+1,childNodes[1],],each=mMax+1) +logB0[,childNodes[2],])
     // adding contribution of all t=n-s>0 values, for s fixed.
     for (s in 0:(mMax-1)){ // s=mMax contributes to n=mMax only, with t=0: already covered
    tMax = mMax-s
    lik[s +1:tMax +1,] = lik[s +1:tMax +1,] +
      exp( lchoose(1:tMax+s,s) + s*logDedge1 + logB0[1:tMax +1,childNodes[1],] + 
           t(logB[s+1,node2sib[childNodes[2]],,1:tMax]))
   }
   // second part:  * (1-Dx)^(-n)
   loglik[,node-nLeaf,] = log(lik) - (0:mMax)*log(1-doomedNode[node])
  }

  /* calculate doomedEdge and transition probs on parent edge.*/
  if (ievent[iedge] == 0){ // BD process
    double gam, psi;
    if (equalBD){
      double tmp = exp(logLam)*edgelength[iedge];
      gam = tmp/(1+tmp);
      psi = gam;
    } else {
      double tmp = exp(lmdiff *edgelength[iedge]);
      psi = (tmp-1)/(tmp-mlratio);
      double gam = psi * mlratio;
    }
    doomedEdge[iedge] = gam + (1-gam)*(1-psi)*doomedNode[iedge]/(1-psi*doomedNode[iedge]);
    double psiprime = 1 - (1-psi)/(1-psi*doomedNode[iedge]);
    double G1 = (1-doomedEdge[iedge])*(1-psiprime);
    // w*(j|i) = G1 w*(j-1|i-1) + psi' w*(j-1|i), for j>=i > 0
    int ie = mM2*iedge;
    for (int js=1; j=mM, jm1=0;  j<mM2;   jm1=j, j+=mM, js++;)
      for (int i=1; i<=js; i++;)
	transitionProb[i+ j+ ie] = G1 * transitionProb[i-1 + jm1 + ie] +
	                       psiprime * transitionProb[i + jm1 + ie];
  }
  if (ievent[iedge] > 0){ // WGD or WGT
    double p1 = p1g[ievent[iedge]-1]; // probs to retain 1, 2 or 3 copies
    double p2 = p2g[ievent[iedge]-1];
    double p3 = p3g[ievent[iedge]-1]; // this is 0 if WGT
    double doo1 = doomedNode[iedge];
    double doo2 = doo1*doo1;
    double doo3 = doo2*doo1;
    doomedEdge[iedge] = p1 * doo1 + p2 * doo2 + p3 * doo3;
    double G1 = (p1 + 2*p2*doo1  +3*p3*doo2) * (1-doo1);
    double G2 = (p2g + 3*p3*doo1) * (1-doo1) * (1-doo1);
    double G3 = p3 * (1-doo1) * (1-doo1) * (1-doo1);
    int ie = mM2*iedge;
    transitionProb[1+ mM*1 + ie]=G1; // w*(1|1)
    transitionProb[1+ mM*2 + ie]=G2; // w*(2|1)
    transitionProb[1+ mM*3 + ie]=G3; // w*(3|1)
    // w*(j|i) = G1 w*(j-1|i-1) + G2 w*(j-2|i-1)      for i>=2     j in [i,2i], or
    // w*(j|i) = G1 w*(j-1|i-1) + G2 w*(j-2|i-1) + G3 w*(j-3|i-1), j in [i,3i]
    for (int i=2; i<mM; i++;)
      for (int j=i*mM, jm1=j-mM, jm2=jm1-mM, jm3=jm2-mM;
	   j<mM2 && ( j<=2*i*mM || (p3>0.0 && j<=3*i*mM) ) ; 
	   jm3=jm2, jm2=jm1, jm1=j, j+=mM;){
	transitionProb[i+ j+ ie] = G1* transitionProb[i-1 + jm1 + ie] +
	                           G2* transitionProb[i-1 + jm2 + ie];
	if (jm3>=0 && p3>0.0)
	  transitionProb[i+j+ie] += G3*transitionProb[i-1 + jm3 + ie];
      }
  }

  /* calculate Bs along parent edge, unless root edge */
  if (ievent[iedge] >= 0){ // non-root
    /* start with logB0: B values for t=0 */
    if (di < *nLeaf){ // leaf edge: B0 = transitionProb(s,data at the tip)
    // fixit: leafName = as.character(phyloMat$Species[edge])
      for (int f=0; f<nF; f++){
        int idata = f + nF * di;
        if (geneCountData[idata]<0){ // missing gene count. B0[i] = (1-doomedEdge)^i
          // fixit to use original NA values
          double tmp=0.0;
          double l1mde = log(1.0-doomedEdge[iedge]);
          for (int is=0, i=mM*di+mMnN*f; is<mM; is++, i++, tmp+=l1mde;)
            logB0[i] = tmp;
        } else {
          for (int i=0 + mM*di + mMnN*f, it= 0+ mM*geneCountData[idata] + mM2*di;
               i<mM; i++, it++;) // B0[i]=transitionProb[i,count]
            logB0[i] = log( transitionProb[it] );
        }
      }
    } else { // internal edge: B0 = transitionProb * lik at descendant node (matrices)
             // i.e. B0[i] = sum transitionProb[i,j] * lik[j], for that node and family.
      for (int f=0; f<nF; f++)
        for (int is=0, i=mM*di+mMnN*f; is<mM; is++, i++){
          for (int js=0, j=mMnN*di, jll=mM*(di- *nLeaf) +mMniN*f;
               js<mM;    js++, j+=mM, jll++;)
            logB0[i] += transitionProb[is+j] * loglik[jll];
          logB0[i] = log(logB0[i]);
        }
    }
    /* logB: B values for t>0, if node is a second sibling */
    if (scdsib[iedge]){
      double lde = log(doomedEdge[iedge]) // log of Gx(0)
      // initialize with case s=mMax: B_t,mMax = Gx(0)^t * B_0,mMax
      // warning: t has index t-1 because logB does not have case t=0.



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

  // fixit: collect contribution to parent node
  }
 return(list(loglikRoot=loglik[,root-nLeaf,],doomedRoot=doomedNode[root],
             doomedRootLeft =doomedEdge[childEdges[1]],doomedRootRight=doomedEdge[childEdges[2]]))
}

# .logsumexp <- function(x,y){ # gets log( exp(x) + exp(y) ) carefully
#  tmp = pmax(x,y); return(tmp + log(1+exp(pmin(x,y)-tmp))) }
