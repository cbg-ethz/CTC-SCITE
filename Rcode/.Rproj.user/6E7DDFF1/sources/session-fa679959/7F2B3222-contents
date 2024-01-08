library(readr)
library(Rcpp)

sourceCpp('mutations_placement.cpp')
find_most_recent_common_ancestor <- function(leaf1, leaf2, ################ Untested, can remove many of the parameters
                                             tree_parent_vector_format,
                                             nCells, nMutations, nSamples,
                                             ancestorMatrix, alleleCount,
                                             leafClusterId,mutatedReadCounts,
                                             totalReadCounts,
                                             dropoutRate, seqErr,
                                             wbcStatus
                                             ){
  ##Trace back the lineage of the tree for one leaf.
  ##Then trace back the lineage of the tree for the other leaf and for every node
  ##whether is lies in the lineage of the first leaf.
  ##The first node that does is the most recent common ancestor node.
  ##Concatenating these two will form the shortest path through the tree.
  ## If there is a mutation on the tree, then this means that the cells are
  ##split by the tree, if there is none, then they aren't.
  
  ##Note that the nodes of the tree are encoded from 0 to the number of nodes -1
  lineage1 <- vector()
  while(lineage1[length(lineage1)]!= length(tree_parent_vector_format))  
    lineage1 <- c(lineage1, tree_parent_vector_format[lineage1[length(lineage1)]])
  }
  lineage2 <- vector()
  while(TRUE)  {
    lineage2 <- c(lineage2, tree_parent_vector_format[lineage2[length(lineage2)]]) 
    if (lineage2[length(lineage2] %in% lineage1) break
  }
  MRCA <- lineage2[length(lineage2)]
  return(c(lineage1,lineage2, MRCA))

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
for(i in 1:nClusters) ClusterID <- c(ClusterID, rep.int(i,description$CellCount[i]))

##Pull apart the count file into counts for mutated read and total counts respectively
mutatedReadCounts <- matrix(0,nrow = nMutations, ncol = 0)
for (j in 1:nClusters){
    print(j)
    mutatedReadCounts <- cbind(mutatedReadCounts,counts[,4+2*j])
}

wildtypeReadCounts <- matrix(0,nrow = nMutations, ncol = 0)
for (j in 1:nClusters){
  print(j)
  wildtypeReadCounts <- cbind(wildtypeReadCounts,counts[,4+2*j-1])
}

totalReadCounts <- mutatedReadCounts + wildtypeReadCounts


##wbc status indicates which of the cells is a white blood cells and which one isn't.
##So far, the cells are arbitrary, and I will assign the fist cells from a cluster to be WBCs.
wbcStatus <- rep(0, nCells)

for(i in 1:nClusters){
  j <- 1
  while(j>0 & j<description$WBCs[i]+1){ #Iterating over the number of White blood cells of a cluster
    wbcStatus[which(ClusterID == i)[1] + j-1] <- 1 # and identifying the first cell
                                              # that belongs to a cluster and counting from then on
  j<- j+1
  }
}



#for (i in nrow(postSampling)){
#  tree <- postSampling$Tree[i]
tree <- postSampling$Tree[1] ##Debugging
  treeParentVectorFormat <- as.numeric(unlist(strsplit(tree, " ")))
  dropoutRate <- postSampling$DropoutRate[i]
  seqErrRate <- postSampling$SequencingErrorRate[i]
  
  ### Now I need to compute the best mutation placement on the tree. This is done
  ##using the scoreTree C++ function (taken from CTC_treeScoring.cpp).
  
  find_most_recent_common_ancestor(1,3, )
  
  ancestorMatrix <- parentVector2ancMatrix(treeParentVectorFormat,
                                            length(treeParentVectorFormat)) #############Something wrong here still
  
  best_mutation_placement <- getMutationPlacement (nCells, nMutations, nSamples, ############# Untested
                                                   ancestorMatrix, alleleCount,
                                                   leafClusterId,mutatedReadCounts,
                                                   totalReadCounts,
                                                   dropoutRate, seqErr, 1,
                                                   wbcStatus 
  
#}


dropoutRate
