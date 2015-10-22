.getMatAndMatDoomed <-
function(logLamlogMu, nLeaf, nFamily, phyloMat, geneCountData, nPos, wgdTab) {
	# phyloMat, wgdTab: read processInput.R for explanation 	
	
	nNode = max(phyloMat$Child)	
	root = nLeaf +1
	
	Mat = array (0, c(nPos, nNode-nLeaf, nFamily))
	# Mat is a 3D matrix of size nPos x (nNode-nLeaf) x nFamily
	# Mat has nNode-nLeaf columns, each column j corresponds to an internal node j+nLeaf
	# each entry Mat[i,j,k] is the prob of the data below node (j+nLeaf) in family k given 
	#   that the node was started with i genes
	
	MatDoomed = array (0, c(nNode,3))
	# MatDoomed is a 2D matrix with nNode rows and 3 columns
	# each column  corresponds to the probability to be doommed when a lineage starts at
	# the Node is considered
	# column 1 probability to be doomed 
	# column 2 probability to be doomed only on the left
	# column 3 probability to be doomed only on the right
	
	
	maxInteNoSingNode = nNode - 2*(nrow(wgdTab))
	 # this is the largest internal-non-singleton node
	 # we calculate lklihd from the maxInteNoSingNode to the root

	#uses the fact that speciation nodes are given lower numbers 
	#than singleton nodes when the tree is read in by the phyext function.

	
	for (parent in maxInteNoSingNode : root) {
		parentIndex = which(phyloMat$Parent == parent)
		 # since parent is in inte-NoSing-Node, length(parentIndex) is 2, i.e. 2 children
		
		child1 = phyloMat$Child[ parentIndex[1]]
		child2 = phyloMat$Child[ parentIndex[2]]
				
		branchNode1 = c(parent, child1)
		branchNode2 = c(parent, child2)
		 # branchNode is a vector of nodes on 1 branch, including singleton nodes at the middle of branch
		 # ex: if no singleton nodes on the branch, length(branchNode) = 2
		 #     if 4 singletons (=2 WGD) on the branch, length(branchNode) = 6
		 
		
		while(child1 > maxInteNoSingNode){
			child1 = phyloMat$Child[phyloMat$Parent==child1]
			branchNode1 = c(branchNode1, child1)
		}
		
		while(child2 > maxInteNoSingNode){
			child2 = phyloMat$Child[phyloMat$Parent==child2]
			branchNode2 = c(branchNode2, child2)
		}
		
		a = .calcProbOneBranch(branchNode1, logLamlogMu, nPos, nFamily, nLeaf, geneCountData, wgdTab, phyloMat, MatDoomed, MDcol=2, Mat)

		prob1 = a$prob; MatDoomed = a$MatDoomed; Mat=a$Mat;

		b = .calcProbOneBranch(branchNode2, logLamlogMu, nPos, nFamily, nLeaf, geneCountData, wgdTab, phyloMat, MatDoomed, MDcol=3, Mat)
		
		prob2 = b$prob; MatDoomed = b$MatDoomed; Mat=b$Mat;
				
		# calcProbOneBranch = function(branchNode, logLamlogMu, nPos, nFamily, nLeaf, 
		#  geneCountData, wgdTab, phyloMat, MatDoomed, MDcol, Mat){
		
		Mat[, parent-nLeaf, ] = prob1 * prob2
		MatDoomed[parent,1]= MatDoomed[parent,2]*MatDoomed[parent,3]		
	}
	return(list(first=Mat,second=MatDoomed))
}
