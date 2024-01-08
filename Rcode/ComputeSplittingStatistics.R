#install.packages("readr")
#install.packages("Rcpp")
#install.packages("dplyr")
library(readr)
library(dplyr)
library(Rcpp)

sourceCpp('mutations_placement.cpp')
find_most_recent_common_ancestor <- function(treeParentVectorFormat, leaf1, leaf2){
  ##Trace back the lineage of the tree for one leaf.
  ##Then trace back the lineage of the tree for the other leaf and for every node
  ##whether is lies in the lineage of the first leaf.
  ##The first node that does is the most recent common ancestor node.
  ##Concatenating these two will form the shortest path through the tree.
  ## If there is a mutation on the tree, then this means that the cells are
  ##split by the tree, if there is none, then they aren't.
  
  ##Note that the nodes and leaves of the tree are encoded from 0 to the number of nodes minus 1
  ## Therefore, I add 1 to the indices to be compatible with R indication starting at 1
  lineage1 <- leaf1
  repeat {  
    #print(treeParentVectorFormat[lineage1[length(lineage1)] + 1])
    lineage1 <- c(lineage1, treeParentVectorFormat[lineage1[length(lineage1)] + 1])
    if(lineage1[length(lineage1)] == length(treeParentVectorFormat)) break
  }
  lineage2 <- leaf2
  nextParent <- treeParentVectorFormat[leaf2 + 1]
  while(!(nextParent %in% lineage1))  {
    lineage2 <- c(lineage2, nextParent)
    nextParent <- treeParentVectorFormat[nextParent + 1]
    #print(nextParent)
    #print(!(nextParent %in% lineage1))
    
  }
  MRCA <- nextParent
  return(list(lineage1,lineage2, MRCA))
}


posteriorSamplingFile <-  "../../input_folder/Br7/Br7_1M_2_seed13543_postSampling.tsv"
countFile <- "../../input_folder/Br7/Br7.txt"
descriptionFile <- "../../input_folder/Br7/Br7_samples_nodeDescription.tsv"

postSampling <- read_delim(posteriorSamplingFile,
                           delim = "\t", col_names = c("LogScore", "SequencingErrorRate","DropoutRate", "LogTau", "Tree"))

counts <- read_delim(countFile,
                     delim = "\t", col_names = FALSE)
description <- read_delim(descriptionFile,
                          delim = "\t", col_names = c("Cluster", "CellCount", "TCs", "WBCs", "Description"))
nCells <- sum(description$CellCount)
nClusters <- nrow(description)  
nMutations <- nrow(counts)
alleleCount <- description$CellCount*2
##The Cluster-ID maps cells identities to the cell-clusters they belong to.
##The i-th entry having value x means that 
## Cells i is in cluster description$Cluster[j]
ClusterID <- vector()
for(i in 1:nClusters) ClusterID <- c(ClusterID, rep.int(i-1,description$CellCount[i]))
## Note that Cpp counts arrays from zero, so the cluster IDs are counted likewise
## in order to be compatible with Cpp code.

##Pull apart the count file into counts for mutated read and total counts respectively
mutatedReadCounts <- matrix(0,nrow = nMutations, ncol = 0)
for (j in 1:nClusters){
    mutatedReadCounts <- cbind(mutatedReadCounts,counts[,4+2*j])
}

wildtypeReadCounts <- matrix(0,nrow = nMutations, ncol = 0)
for (j in 1:nClusters){
  wildtypeReadCounts <- cbind(wildtypeReadCounts,counts[,4+2*j-1])
}


totalReadCounts <- mutatedReadCounts + wildtypeReadCounts


mutatedReadCounts <- mutatedReadCounts %>% t() %>% as.data.frame() %>% as.list()
wildtypeReadCounts <- wildtypeReadCounts %>% t() %>% as.data.frame() %>% as.list()
totalReadCounts <- totalReadCounts %>% t() %>% as.data.frame() %>% as.list()



##wbc status indicates which of the cells is a white blood cells and which one isn't.
##So far, the cells are arbitrary, and I will assign the fist cells from a cluster to be WBCs.
wbcStatus <- rep(0, nCells)

for(i in 1:nClusters){
  j <- 1
  while(j <= description$WBCs[i]){ #Iterating over the number of White blood cells of a cluster
    wbcStatus[which(ClusterID == i-1)[1] + j-1] <- 1 # and identifying the first cell
                                              # that belongs to a cluster and counting from then on
  ## Note: The cluster IDs are counted from zero!  
    j<- j+1
  }
}



#for (i in nrow(postSampling)){
#  tree <- postSampling$Tree[i]
tree <- postSampling$Tree[3200] ##Debugging
  treeParentVectorFormat <- as.numeric(unlist(strsplit(tree, " ")))
  dropoutRate <- postSampling$DropoutRate[i]
  seqErrRate <- postSampling$SequencingErrorRate[i]
  
  ### Now I need to compute the best mutation placement on the tree. This is done
  ##using the scoreTree C++ function (taken from CTC_treeScoring.cpp).
  
  
  ancestorMatrix <- parentVector2ancMatrix(treeParentVectorFormat,
                                            length(treeParentVectorFormat)) #%>%
#    matrix(ncol = length(treeParentVectorFormat), byrow = TRUE)
  
  sourceCpp('mutations_placement.cpp') 
  

  bestMutationPlacement <- getMutationPlacement (nCells, nMutations, nClusters,
                                                   ancestorMatrix, alleleCount,
                                                   ClusterID,mutatedReadCounts,
                                                   totalReadCounts,
                                                   dropoutRate, seqErrRate, 1,
                                                   wbcStatus)

  pairwiseGenealogy <- find_most_recent_common_ancestor(treeParentVectorFormat, 0,2)
  pairwiseGenealogy
  
  positionOfMRCA <- which(pairwiseGenealogy[[1]] == pairwiseGenealogy[[3]])
  firstLeafToMRCA <- pairwiseGenealogy[[1]][1:(positionOfMRCA)]
  
  secondLeafToMRCA <- pairwiseGenealogy[[2]]
  
  #print(firstLeafToMRCA)
  #print(secondLeafToMRCA)
  pathBetweenLeaves <- c(firstLeafToMRCA,rev(secondLeafToMRCA))
  #print(pathBetweenLeaves)
  
  ## Now count how many of the mutations lie in the shortest path between the leaves.
  ##Need to exclude the MRCA for this

  sum(bestMutationPlacement %in% pathBetweenLeaves[pathBetweenLeaves != firstLeafToMRCA[positionOfMRCA]])
  
#}

