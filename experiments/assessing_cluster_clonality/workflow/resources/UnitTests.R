############
# Test functions
############


test_compute_pairwise_distance_of_leaves1 <- function() {
  treeParentVectorFormat <- c(6, 8, 1, 1, 3, 3, 8, 6)
  bestMutationPlacement <- c(5, 8, 6, 0)

  outcome <- c(0, 0)


  pairwiseGenealogy <- findMostRecentCommonAncestor(treeParentVectorFormat, 5, 5)


  dist <- computePairwiseDistanceOfLeaves2(
    treeParentVectorFormat, 5, 5,
    bestMutationPlacement, pairwiseGenealogy
  )

  PASS <- TRUE

  if (sum(dist != outcome)) {
    PASS <- FALSE
  }
  return(PASS)
}

test_compute_pairwise_distance_of_leaves2 <- function() {
  treeParentVectorFormat <- c(6, 8, 1, 1, 3, 3, 8, 6)
  bestMutationPlacement <- c(5, 8, 6, 0)

  outcome <- c(0, 0)

  pairwiseGenealogy <- findMostRecentCommonAncestor(treeParentVectorFormat, 2, 4)

  dist <- computePairwiseDistanceOfLeaves2(
    treeParentVectorFormat, 2, 4,
    bestMutationPlacement, pairwiseGenealogy
  )

  PASS <- TRUE

  if (sum(dist != outcome)) {
    PASS <- FALSE
  }
  return(PASS)
}



test_compute_pairwise_distance_of_leaves3 <- function() {
  treeParentVectorFormat <- c(6, 8, 1, 1, 3, 3, 8, 6)
  bestMutationPlacement <- c(5, 8, 6, 0)
  outcome <- c(3, 1)

  pairwiseGenealogy <- findMostRecentCommonAncestor(treeParentVectorFormat, 0, 5)

  dist <- computePairwiseDistanceOfLeaves2(
    treeParentVectorFormat, 0, 5,
    bestMutationPlacement, pairwiseGenealogy
  )

  PASS <- TRUE

  if (sum(dist != outcome)) {
    PASS <- FALSE
  }
  return(PASS)
}





test_find_MRCA1 <- function() {
  treeParentVectorFormat <- c(6, 8, 1, 1, 3, 3, 8, 6)


  res <- findMostRecentCommonAncestor(treeParentVectorFormat, 5, 5)

  outcome <- list(c(5, 3, 1, 8), c(), 5, 5)
  PASS <- TRUE

  for (i in c(1, 3, 4)) {
    if (sum(outcome[[i]] != res[[i]]) > 0) {
      PASS <- FALSE
      break
    }
  }

  if (length(res[[2]]) > 0) {
    PASS <- FALSE
  }


  return(PASS)
}

test_find_MRCA2 <- function() {
  treeParentVectorFormat <- c(6, 8, 1, 1, 3, 3, 8, 6)
  res <- findMostRecentCommonAncestor(treeParentVectorFormat, 2, 4)

  outcome <- list(c(2, 1, 8), c(4, 3), 1, c(2, 1, 3, 4))
  PASS <- TRUE

  for (i in 1:4) {
    if (sum(outcome[[i]] != res[[i]]) > 0) {
      PASS <- FALSE
      break
    }
  }



  return(PASS)
}

test_find_MRCA3 <- function() {
  treeParentVectorFormat <- c(6, 8, 1, 1, 3, 3, 8, 6)
  res <- findMostRecentCommonAncestor(treeParentVectorFormat, 1, 6)

  outcome <- list(c(1, 8), 6, 8, c(1, 8, 6), c(1, 8, 6))
  PASS <- TRUE

  for (i in 1:3) {
    if (sum(outcome[[i]] != res[[i]]) > 0) {
      PASS <- FALSE
      break
    }
  }



  return(PASS)
}





test_computePairwiseDistanceOfLeavesGivenTree <- function() {
  source("functions.R")

  inputFolder <- "../../input_folder"
  treeName <- "Br7"

  input <- load_data(inputFolder, treeName)

  postSampling <- input$postSampling

  nCells <- 24
  nMutations <- 10
  nClusters <- 11

  alleleCount <- c(8, 6, 2, 6, 2, 2, 2, 2, 8, 4, 6)
  ClusterID <- c(0, 0, 0, 0, 1, 1, 1, 2, 3, 3, 3, 4, 5, 6, 7, 8, 8, 8, 8, 9, 9, 10, 10, 10)
  mutatedReadCounts <- list(
    c(0, 4, 0, 0, 0, 0, 4, 0, 11, 0, 0),
    c(4, 0, 1, 3, 0, 0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 25, 0, 0, 0, 0, 4, 0, 0),
    c(3, 0, 0, 3, 0, 0, 0, 0, 1, 0, 0),
    c(8, 1, 16, 4, 73, 0, 0, 0, 5, 0, 0),
    c(0, 0, 9, 0, 0, 0, 26, 0, 8, 0, 0),
    c(0, 0, 2, 0, 22, 0, 0, 0, 8, 0, 0),
    c(12, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0),
    c(7, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0),
    c(16, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0)
  )

  totalReadCounts <- list(
    c(31, 238, 0, 234, 0, 0, 20, 0, 147, 0, 245),
    c(16, 31, 8, 68, 0, 0, 3, 0, 34, 0, 0),
    c(5, 5, 64, 114, 0, 0, 128, 0, 51, 0, 0),
    c(13, 16, 13, 26, 10, 0, 15, 0, 14, 3, 14),
    c(120, 98, 102, 503, 181, 0, 0, 0, 50, 0, 156),
    c(22, 6, 47, 0, 0, 0, 64, 0, 23, 14, 5),
    c(0, 63, 62, 0, 45, 0, 171, 0, 99, 0, 93),
    c(24, 12, 0, 6, 2, 0, 0, 0, 7, 0, 0),
    c(14, 0, 0, 0, 0, 0, 0, 0, 49, 0, 24),
    c(32, 42, 0, 282, 0, 0, 0, 0, 119, 0, 19)
  )

  seqErrRate <- 0.00225
  dropoutRate <- 0.301
  wbcStatus <- c(1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1)


  res <- computePairwiseDistanceOfLeavesGivenTree(postSampling$"1", 7, 1,
    nCells, nMutations, nClusters, alleleCount,
    ClusterID, mutatedReadCounts, totalReadCounts, wbcStatus,
    nSamplingEvents = 10
  )
  outcome <- c(4, 7, 3, 0, 0, 0, 0, 0, 2, 0)

  PASS <- TRUE

  if (sum(res != outcome) > 0) PASS <- FALSE

  return(PASS)
}



test_transposeMatrix <- function() {
  matrix <- list(c(1, 2), c(3, 4), c(5, 6))

  res <- transposeMatrix(matrix, 3, 2)

  expected <- list(c(1, 3, 5), c(2, 4, 6))
  PASS <- TRUE

  if (sum(expected[[1]] != res[[1]]) > 0) PASS <- FALSE
  if (sum(expected[[2]] != res[[2]]) > 0) PASS <- FALSE


  return(PASS)
}




test_mutation_distribution <- function() {
  nCells <- 24
  nMutations <- 10
  nClusters <- 11

  alleleCount <- c(8, 6, 2, 6, 2, 2, 2, 2, 8, 4, 6)
  ClusterID <- c(0, 0, 0, 0, 1, 1, 1, 2, 3, 3, 3, 4, 5, 6, 7, 8, 8, 8, 8, 9, 9, 10, 10, 10)
  mutatedReadCounts <- list(
    c(0, 4, 0, 0, 0, 0, 4, 0, 11, 0, 0),
    c(4, 0, 1, 3, 0, 0, 0, 0, 1, 0, 0),
    c(0, 0, 0, 25, 0, 0, 0, 0, 4, 0, 0),
    c(3, 0, 0, 3, 0, 0, 0, 0, 1, 0, 0),
    c(8, 1, 16, 4, 73, 0, 0, 0, 5, 0, 0),
    c(0, 0, 9, 0, 0, 0, 26, 0, 8, 0, 0),
    c(0, 0, 2, 0, 22, 0, 0, 0, 8, 0, 0),
    c(12, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0),
    c(7, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0),
    c(16, 0, 0, 0, 0, 0, 0, 0, 11, 0, 0)
  )

  totalReadCounts <- list(
    c(31, 238, 0, 234, 0, 0, 20, 0, 147, 0, 245),
    c(16, 31, 8, 68, 0, 0, 3, 0, 34, 0, 0),
    c(5, 5, 64, 114, 0, 0, 128, 0, 51, 0, 0),
    c(13, 16, 13, 26, 10, 0, 15, 0, 14, 3, 14),
    c(120, 98, 102, 503, 181, 0, 0, 0, 50, 0, 156),
    c(22, 6, 47, 0, 0, 0, 64, 0, 23, 14, 5),
    c(0, 63, 62, 0, 45, 0, 171, 0, 99, 0, 93),
    c(24, 12, 0, 6, 2, 0, 0, 0, 7, 0, 0),
    c(14, 0, 0, 0, 0, 0, 0, 0, 49, 0, 24),
    c(32, 42, 0, 282, 0, 0, 0, 0, 119, 0, 19)
  )

  seqErrRate <- 0.00225
  dropoutRate <- 0.301
  wbcStatus <- c(1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1)

  tree <- c(34, 44, 26, 41, 33, 30, 38, 35, 37, 31, 26, 43, 44, 40, 32, 45, 25, 35, 42, 24, 39, 32, 28, 27, 30, 33, 42, 29, 46, 46, 25, 37, 39, 38, 28, 40, 31, 24, 45, 27, 36, 36, 43, 41, 34, 29)

  ancestorMatrix <- parentVector2ancMatrix(tree, length(tree))


  mutationDistributions <- computeMutationDistribution(
    nCells, nMutations, nClusters, ancestorMatrix,
    alleleCount, ClusterID, mutatedReadCounts,
    totalReadCounts, dropoutRate, seqErrRate, 1,
    wbcStatus
  )

  PASS <- TRUE
  if (length(mutationDistributions) != nMutations) PASS <- FALSE


  numberOfNodesInTree <- 2 * nCells - 1

  for (i in 1:length(mutationDistributions)) {
    if (length(mutationDistributions[[i]]) != numberOfNodesInTree) PASS <- FALSE
  }


  return(list("logMutationDistributions" = mutationDistributions, "PASS" = PASS))
}





test_transposeMatrix <- function() {
  matrix <- list(c(1, 2), c(3, 4), c(5, 6))

  res <- transposeMatrix(matrix, 3, 2)

  expected <- list(c(1, 3, 5), c(2, 4, 6))
  PASS <- TRUE

  if (sum(expected[[1]] != res[[1]]) > 0) PASS <- FALSE
  if (sum(expected[[2]] != res[[2]]) > 0) PASS <- FALSE


  return(PASS)
}


test_sampleMutationPlacements <- function() {
  nSamplings <- 100000
  nMutations <- 10
  nCells <- 24
  logProbs <- test_mutation_distribution()$logMutationDistributions

  mutationSampling <- transposeMatrix(
    sampleMutationsPlacement(
      nSamplings,
      nMutations, logProbs
    ),
    nSamplings, nMutations
  )


  PASS <- TRUE
  for (mutation in 1:length(mutationSampling)) {
    sampledPlacement <- vector()
    for (i in 0:(2 * nCells - 2)) {
      sampledPlacement <- c(sampledPlacement, sum(mutationSampling[[mutation]] == i))
    }
    sampledPlacement <- sampledPlacement / sum(sampledPlacement)
    Probs <- logProbs[[mutation]] %>% exp()
    Probs <- Probs / sum(Probs)
    if (((Probs - sampledPlacement)^2 %>% sum()) > 10e-4) PASS <- FALSE
  }
  return(PASS)
}


test_ComputePerMutationProbabilityOfPolyclonality <- function() {
  treeParentVectorFormat <- c(6, 8, 1, 1, 3, 3, 8, 6)


  pairwiseGenealogy <- list(c(0, 6, 8), c(2, 1), c(8), c(0, 6, 8, 1, 2))
  nCells <- 5
  nMutations <- 2


  PASS <- TRUE

  logMutationPlacementProbs <- list(log(c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0)), log(c(0, 0.5, 0, 0.5, 0, 0, 0, 0, 0)))
  outcome <- 0.5
  res <- ComputePerMutationProbabilityOfPolyclonality(
    pairwiseGenealogy,
    logMutationPlacementProbs,
    nMutations, nCells
  )
  if (res != outcome) PASS <- FALSE


  logMutationPlacementProbs <- list(log(c(0, 0.6, 0.4, 0, 0, 0, 0, 0, 0)), log(c(0, 0.4, 0, 0.6, 0, 0, 0, 0, 0)))
  outcome <- 0
  res <- ComputePerMutationProbabilityOfPolyclonality(
    pairwiseGenealogy,
    logMutationPlacementProbs,
    nMutations, nCells
  )
  if (res != outcome) PASS <- FALSE

  logMutationPlacementProbs <- list(log(c(0.4, 0.6, 00, 0, 0, 0, 0, 0, 0)), log(c(0, 0.4, 0, 0.6, 0, 0, 0, 0, 0)))
  outcome <- 0.4
  res <- ComputePerMutationProbabilityOfPolyclonality(
    pairwiseGenealogy,
    logMutationPlacementProbs,
    nMutations, nCells
  )
  if (res != outcome) PASS <- FALSE

  return(PASS)
}
