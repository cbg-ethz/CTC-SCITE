source("functions.R")

############
#Config
############
inputFolder <- "../../input_folder"
treeName <- "Br7"




############
#Data preprocessing
############
input <- load_data(inputFolder, treeName)

postSampling <- input$postSampling
nClusters <- input$nClusters
ClusterID <- input$clusterID
nCells <- input$nCells  
nMutations <- input$nMutations
nClusters <- input$nClusters
alleleCount <- input$alleleCount
ClustserID <- input$ClusterID
mutatedReadCounts <- input$mutatedReadCounts
totalReadCounts <- input$totalReadCounts
wbcStatus <- input$wbcStatus



############
#Unit testing
############



test_compute_pairwise_distance_of_leaves1()
test_compute_pairwise_distance_of_leaves2()
test_compute_pairwise_distance_of_leaves3()
test_find_MRCA1()
test_find_MRCA2()
test_find_MRCA3()
test_computePairwiseDistanceOfLeavesGivenTree()
test_transposeMatrix()
############
#Main Analysis
############

desired_values <- sample(1:length(postSampling), size = 4000, replace = FALSE) %>% sort()

postSampling <- postSampling[desired_values]

## Go through all clusters and compare all pairs of cells within each cluster with
## each other. Note that the cells from the clusters are adjacent to each other by
## design, so incrementing the index j by 1 makes sense
distance <- vector()
clusterIdentityofdistance <- vector()
system.time(for (c in 1:nClusters){
  cellsInCluster <- which(ClusterID == (c-1))-1 ## Make sure array indication is 
                                                ## compatible with cpp
  cluster_done <- 0
  for(i in cellsInCluster){
    if(cluster_done == 1){
      cluster_done <- 0
      break
    }
    if(wbcStatus[i+1] == 1) next
    j <- cellsInCluster[1]
    while(j < i){
      if(cluster_done == 1){
        break
      }
      if(wbcStatus[j+1] == 1){
        j <- j + 1
        next
      }
      print(paste(paste("Computing genomic distances of leaves:", i, sep = " "), j, sep = " "))
      distance <- c(distance, produce_Distance_Posterior(i,j,postSampling, treeName, nCells,
                                                         nMutations, nClusters,
                                                         alleleCount,ClusterID,
                                                         mutatedReadCounts, totalReadCounts,wbcStatus, nSamplingEvents = 1000))
      clusterIdentityofdistance <- c(clusterIdentityofdistance, c-1)
      j <- j + 1
      cluster_done <- 1
    }
  }
  
})
#########


intraClusterSplitMedianPlot <- ggplot(data.frame(Median_Distance = distance), aes(x = Median_Distance)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black", alpha = 0.7)+ 
  xlab("Median distance between of leaves within the same cluster") + ylab("total count") +
  ggtitle(treeName) +
  labs(subtitle = "Histogram of similarities of cells within cluster",caption = "dashed red line indicates cutoff for oligoclonality") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    plot.subtitle = element_text(size= 20),
    axis.text = element_text(size = 16),
    plot.caption = element_text(size = 14)
  )

plot(intraClusterSplitMedianPlot)

summary(distance)
print(intraClusterSplitMedianPlot)
#####Manually adapt
cutoffForOligoclonality <- 400


intraClusterSplitMedianPlot + geom_vline(xintercept = cutoffForOligoclonality,color = "red", linetype = "dashed", size = 1)



### Now look at each cluster and determine whether at least one pair of cells split
clusterSplits <- vector()

for (c in 1:nClusters){
  cellPairsInCluster <- which(clusterIdentityofdistance == (c-1))
  print(cellPairsInCluster)
  print("Distances:")
  print(distance[cellPairsInCluster])
  
  if(length(cellPairsInCluster) == 0) next
  else if (max(distance[cellPairsInCluster])>cutoffForOligoclonality){
    clusterSplits <- c(clusterSplits,1)
  }
  else {
    clusterSplits <- c(clusterSplits,0)
  }
}
which(ClusterID == 25)

produce_Distance_Posterior(35,36, postSampling, "LM2")

produce_Distance_Posterior(35,36, postSampling, "LM2")


#compute_hamming_distance_distr <- function(leaf1, leaf2, postSampling){
  
#  for (i in nrow(postSampling)){
#    tree <- postSampling$Tree[i]
    #tree <- postSampling$Tree[3200] ##Debugging
    
#  }  
#}


