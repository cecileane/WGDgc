rgenecount <- function(nfam,tre,lambdamu,retention,geomMean=NULL,dirac=NULL,
                       conditioning=c("none")){ # c("none","oneOrMore","twoOrMore","oneInBothClades")){
  # generates gene count data.
  # nfam: number of gene families to simulate, tre: species tree in SIMMAP format
  # lambdamu: vector of size 1 (to assume lambda=mu) or size 2 (lambda,mu)

  # warning: conditioning *not* implemented yet. If desired, the user needs to
  #          filter the simulated data and re-simulate more families to
  #          get the desired number of families.  

  if (length(lambdamu)==1){
    equalBDrates <- TRUE
    if (lambdamu <= 0) stop("lambdamu must be positive")
    lambda <- lambdamu
  } else if (length(lambdamu)==2) {
    if (identical(TRUE,all.equal(lambdamu[1],lambdamu[2],tolerance = 1e-06))){
      equalBDrates <- TRUE
      lambda <- lambdamu
      if (lambda <= 0) stop("lambdamu must be positive")
    } else{
      equalBDrates <- FALSE
      if (any(lambdamu <= 0)) stop ("lambda and mu must be positive")
      lmdiff <- lambdamu[1] - lambdamu[2]  # lambda - mu
      mlratio <- lambdamu[2] / lambdamu[1] # mu/lambda
      #cat("lmdiff=",lmdiff,"mlratio=",mlratio,"\n")
    }
  } else stop("lambdamu must be of length 1 or 2.")

  if( (is.null(geomMean) + is.null(dirac)) != 1 )
    stop("Use exactly one: geomMean or dirac.")

  input <- processInput(tre,equalBDrates=equalBDrates,startingQ=retention)

  nWGD=sum(input$phyloMat$type == "WGD")
  nWGT=sum(input$phyloMat$type == "WGT")
  nWG = nWGD + nWGT
  if (length(retention) != nWG)
    stop("retention needs to have length equal to the number of WG events.")

  conditioning <- match.arg(conditioning)

  alldat <- matrix(NA,nfam,input$nNode) # simulated data, at all nodes, one column per row

  #----- starting simulations: root node --------------#
  node <- input$nLeaf+1
  if (!is.null(dirac)){
    if (dirac != floor(dirac)){
      dirac <- floor(dirac)
      warning("the dirac value was rounded down to be an integer.",immediate.=TRUE)
    }
    if (dirac<1) stop("the dirac value needs to be positive.")
    alldat[,node] <- dirac
  }
  geomProb <- NULL
  if (!is.null(geomMean)){
    if (geomMean<1)
      stop("the mean number of ancestral genes (geomMean) has to be >= 1.")
    geomProb <- 1/geomMean
    alldat[,node] <- rgeom(nfam,prob=geomProb) + 1
  }

  #----- continuing simulations: along all non-root edges
  # reversed post-order: root to tips. Last edge was root edge.
  for (ie in (nrow(input$edgeOrder)-1):1){
    node <- input$edgeOrder$child[ie] # child node
    edge <- input$edgeOrder$edge[ie]
    x0 <- alldat[,input$phyloMat$Parent[edge]] # number of starting gene lineages
    #cat("Starting edge",edge,"\n")

    if (input$edgeOrder$type[ie] == "BD"){
      len = input$phyloMat$Time[edge] # branch length
      # gam = transition prob from 1 to 0.
      # psi = probability parameter for transition from 1 to k, k>=1: (1-gam)(1-psi) psi^{k-1}
      if (equalBDrates){
        gam <- lambda*len/(1+lambda*len)
        psi <- gam
      } else {
        psi <- (exp(lmdiff *len) -1) / (exp(lmdiff *len) -mlratio)
        gam <- psi * mlratio
      }
      #cat("  BD edge, length",len,"gam=",gam,"psi=",psi,"\n")
      x1 <- rbinom(nfam,size=x0,prob=1-gam)   # how many lineages lead to >0 genes
      alldat[x1==0,node] <- 0
      # Warning: if this procedure is not efficient, need to implement
      #          the calculation of transition probabilities, from any m to n.
      #          Using m=1 only below.
      ind <- which(x1>0) # families with 1 or more descendants
      l <- length(ind)
      x1 <- x1[ind]
      endlin <- cumsum(x1)
      z <- rgeom(endlin[l],prob=1-psi) # z+1 = end number of daughter lineages
      startlin <- 1+ c(0, endlin[1:(l-1)])
      for (i in 1:length(x1))
        alldat[ind[i],node] <- sum(z[startlin[i]:endlin[i]]) + x1[i]
    }
    if (input$edgeOrder$type[ie] == "WGD"){
      ievent <- which(input$wgdTab$wgdNode == input$phyloMat$Parent[edge])
      if (length(ievent) != 1) stop("more or less than 1 event matched the WGD node.")
      p1g <- input$wgdTab$retain1[ievent] # prob to retain 1 gene (= only the original copy)
      p2g <- input$wgdTab$retain2[ievent] #                2 genes
      cat("  WGD edge, probs",p1g,"and",p2g,"\n")
      alldat[,node] <- x0 + rbinom(nfam,size=x0,prob=p2g)
    }
    if (input$edgeOrder$type[ie] == "WGT"){
      ievent <- which(input$wgdTab$wgdNode == input$phyloMat$Parent[edge])
      if (length(ievent) != 1) stop("more or less than 1 event matched the WGD node.")
      p1g <- input$wgdTab$retain1[ievent]
      p2g <- input$wgdTab$retain2[ievent]
      p3g <- input$wgdTab$retain3[ievent]
      cat("  WGT edge, probs",p1g,",",p2g,"and",p3g,"\n")
      # first simulate the number of genes (x1) that retain only 1 copy
      # then look at the x0-x1 genes that retain either 2 or 3:
      #           of those, simulate the number x2 that retain only 2.
      # final number: x1 + 2*x2 + 3*(x0-x1-x2) = x0 + x2 + 2*(x0-x1-x2) = x0 - x2 + 2*(x0-x1)
      x1 <- rbinom(nfam,size=x0,   prob=p1g)
      x2 <- rbinom(nfam,size=x0-x1,prob=p2g/(p2g+p3g))
      alldat[,node] <- x0 - x2 + 2*(x0-x1)
    }
  }
  # fixit: implement conditioning
  #if (conditioning=="oneInBothClades") {
  #  .checkClades(input$phyloMat, geneCountData, input$nLeaf)
  colnames(alldat) <- 1:input$nNode
  for (isp in 1:input$nLeaf){
    edge <- which(input$phyloMat$Child == isp)
    colnames(alldat)[isp] <- as.character(input$phyloMat$Species[edge])
  }
  colnames(alldat)[input$nLeaf+1] <- "root"
  return(alldat)
}
