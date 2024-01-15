
############
#Function Definitions
############

library(readr)
library(dplyr)
library(Rcpp)
library(purrr)
library(ggplot2)

sourceCpp('mutations_placement.cpp')


#Input: a tree in parent vector format, meaning that the i'th entry of the vector
# is te parent node of the entry i. Nodes are counted from zero and the root is 
# length(Tree)
#Output: A list with three entries:
#  - the first entry is a vector of nodes tracing back leaf 1 to the root
#  - the second entry is a vector of nodes tracing back leaf2 to the descendant
#    of the MRCA in the lineage
#  - the thirs entry is the MRCA

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

# Take a tree and a pair of mutations and do the following:
compute_pairwise_distance_of_leaves <- function(treeData, leaf1, leaf2, nCells,
                                                nMutations, nClusters,
                                                alleleCount,ClusterID,
                                                mutatedReadCounts, totalReadCounts,wbcStatus){
  tree <- treeData$Tree
  treeParentVectorFormat <- as.numeric(unlist(strsplit(tree, " ")))
  dropoutRate <- treeData$DropoutRate
  seqErrRate <- treeData$SequencingErrorRate
  
  ### Now I need to compute the best mutation placement on the tree. This is done
  ##using the scoreTree C++ function (taken from CTC_treeScoring.cpp).
  
  #print("Preprocess tree")
  ancestorMatrix <- parentVector2ancMatrix(treeParentVectorFormat,
                                           length(treeParentVectorFormat))
  
  
  #print("Find best Mutation placement")
  bestMutationPlacement <- getMutationPlacement (nCells, nMutations, nClusters,
                                                 ancestorMatrix, alleleCount,
                                                 ClusterID,mutatedReadCounts,
                                                 totalReadCounts,
                                                 dropoutRate, seqErrRate, 1,
                                                 wbcStatus)## This crashes when executed on posterior sampling of Br61
  #print("Finding most recent common ancestor")
  pairwiseGenealogy <- findMostRecentCommonAncestor(treeParentVectorFormat, leaf1,leaf2)
  
  positionOfMRCA <- which(pairwiseGenealogy[[1]] == pairwiseGenealogy[[3]])
  firstLeafToMRCA <- pairwiseGenealogy[[1]][1:(positionOfMRCA)]
  
  secondLeafToMRCA <- pairwiseGenealogy[[2]]
  
  #print(firstLeafToMRCA)
  #print(secondLeafToMRCA)
  pathBetweenLeaves <- c(firstLeafToMRCA,rev(secondLeafToMRCA))
  #print(pathBetweenLeaves)
  
  ## Now count an output how many of the mutations lie in the shortest path between the leaves.
  ##This equals the Hamming distance between the inferred Genotypes of two leaves
  ##Need to exclude the MRCA for this
  
  return(sum(bestMutationPlacement %in% pathBetweenLeaves[pathBetweenLeaves != firstLeafToMRCA[positionOfMRCA]]))
}

###The function computes the distribution of evolutionary distances of two
##specified leaves as sampled from the posterior distribution of trees.
## @Input: leaf1: integer-valued index of first leaf
## @Input: postSampling: loaded tibble containing the posterior Sampling
## @Input: treeName: 
## @Output: ....



produce_Distance_Posterior <- function(leaf1, leaf2,postSampling, treeName,nCells,
                                       nMutations, nClusters,
                                       alleleCount,ClusterID,
                                       mutatedReadCounts, totalReadCounts,wbcStatus){
  
  ## For each row in the posterior Sampling file, the distance of two leaves is computed

  
  dist_histogram <- sapply(split(postSampling, seq(nrow(postSampling))),
                           FUN = compute_pairwise_distance_of_leaves, leaf1,leaf2,
                           nCells, nMutations,nClusters, alleleCount,
                           ClusterID, mutatedReadCounts, totalReadCounts, wbcStatus)
  median_dist <- median(dist_histogram)
  
  
  plot(
    ggplot(data.frame(HammingDistance = dist_histogram), aes(x = HammingDistance)) +
      geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7)+ 
      xlab(sprintf("genetic distance between leaf %d and leaf %d", leaf1, leaf2)) + ylab("total count") +
      ggtitle(paste("Histogram of genetic distances of clusters cells in", treeName)) +
      geom_vline(xintercept = median_dist,color = "red", linetype = "dashed", size = 1) +
      labs(subtitle = "As sampled from the posterior distribution",caption = "median indicated by dashed red line") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 24, face = "bold"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.subtitle = element_text(size= 20),
        axis.text = element_text(size = 16) 
      )
  )
  
  
  return(median(dist_histogram))
}

# Loads all necessary data for the CTC-project.
# Specifically it return a named list as follows:
# postSampling: Loads the posterior sampling tsv as a tibble with the 
#               following columns:
#               the (unnormalised) LogScore, estimated sequencing error rate,
#               the estimated dropout rate, logTau and the Tree in parent vector 
#               format meaning that the i'th entry of the vector
#               is te parent node of the entry i. Nodes are counted from zero and the root is 
#               length(Tree)
# nClusters: The total number of CTC-clusters
# ClusterID: The Cluster-ID maps cells identities to the cell-clusters they belong to.
#            The i-th entry having value x means that 
#            Cells i is in the j-th cluster (7th row in the description file)
# nCells:     total number of Cells in the experiment
# nMutations: Total number of considered mutations in the dataset
# nClusters: Total number of CTC-clusters in the experiment
# alleleCount: Total number of cells per Cluster*2
# ClusterID: Number 
# mutatedReadCounts: A vector containing the number of of mutated reads for each
#                    cluster
# total Read Counts: vector containingtotal read count for each cluster
# wbcStatus: A vector that has value 1 if cell i is a white blood cells and 
#            0 else


load_data <- function(inputFolder, treeName){
  ##Define paths
  
  posteriorSamplingFile <-  sprintf("%s/%s/%s_postSampling.tsv", inputFolder, treeName,treeName)
  
  countFile <- sprintf("%s/%s/%s.txt", inputFolder, treeName,treeName)
  descriptionFile <- sprintf("%s/%s/%s_samples_nodeDescription.tsv", inputFolder, treeName, treeName)
  
  
  ##Load data
  
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
  
  

  
  
  
  return(list("postSampling" = postSampling, "nClusters" = nClusters,
              "clusterID" = ClusterID, "nCells" = nCells,
              "nMutations" = nMutations, "nClusters" = nClusters,
              "alleleCount" = alleleCount, "ClusterID" = ClusterID,
              "mutatedReadCounts" = mutatedReadCounts,
              "totalReadCounts" = totalReadCounts, "wbcStatus" = wbcStatus))
}


