############
# Function Definitions
############

library(Rcpp)
library(tidyverse)



sourceCpp("../../workflow/resources/mutations_placement.cpp")



#' Takes a list of mutations and outputs which one of these is a driver.
#'
#'
#' @param mutations a names vector containing chromosomes in the format "chrN"
#'                  in the first column and an integer chromosomal position on
#'                  the second column
#' @param annotations an annotation data frame. Must contain the columns
#'  - 'CGI-Oncogenic Summary': entry can be 'driver (oncodriveMUT)' or
#'    somerthing else
#'  - 'CGI-Oncogenic Prediction': entry can be 'oncogenic (predicted)' or
#'    something else
#'  - 'CGI-External oncogenic annotation'
#'
#' @return a boolean vector with as many entries as there are rows in mutations
#' @export
#'
#' @examples
IsDriver <- function(mutations, annotations) {
  annotated_mutations <- annotations %>%
    filter(annotations$"#CHROM" == as.character(mutations[1]) &
      annotations$POS == as.numeric(mutations[2]))

  check <- annotated_mutations %>%
    select(c(
      "CGI-Oncogenic Summary", "CGI-Oncogenic Prediction",
      "CGI-External oncogenic annotation"
    )) %in%
    c("oncogenic (predicted)", "driver (oncodriveMUT)") %>%
    sum()

  driver <- FALSE
  if (check > 0) {
    driver <- TRUE
  }
  return(driver)
}




#' computes the distribution of evolutionary distances of two
#' specified leaves as sampled from the posterior distribution of trees.
#'
#' @param leaf1 integer-valued index of first leaf
#' @param leaf2 integer-valued index of second leaf
#' @param postSampling loaded list of tibbles tibble containing the posterior
#'                     Sampling
#' @param treeName character string: Name of the tree for the output plot
#' @param nCells integer-valued total number of cells in the dataset
#' @param nMutations integer-valued total number of mutations in the dataset
#' @param nClusters integer-valued total number of clusters in the dataset
#' @param alleleCount integer vector of numbers of alleles per clusters
#' @param ClusterID integer vector of cluster IDs
#' @param mutatedReadCounts list of integer-valued vectors indicating the number
#' of mutated read per mutation (list index) and sample (vector index)
#' @param totalReadCounts list of integer-valued vectors indicating the total
#' number of reads per mutation (list index) and sample (vector index)
#' @param wbcStatus boolean vector of length nCells indicating for each cell if
#' it is a white blood cell (TRUE) or not (FALSE)
#'
#' @return splittingFraction: The fraction of sampling events for which the pair
#'         of cells
#' shows branching evolution
#'
#' @export
#'
#' @examples
produce_Distance_Posterior <- function(leaf1, leaf2, postSampling, treeName,
                                       nCells, nMutations, nClusters,
                                       alleleCount, ClusterID,
                                       mutatedReadCounts, totalReadCounts,
                                       wbcStatus, nSamplingEvents = 20,
                                       clusterName = "") {
  ## For each row in the posterior Sampling file, the distance of two leaves is
  ## computed

  print("Computing the posterior distribution")

  distance_statistics <- parallel::mclapply(postSampling,
    FUN = computePairwiseDistanceOfLeavesGivenTree, leaf1, leaf2,
    nCells, nMutations, nClusters, alleleCount,
    ClusterID, mutatedReadCounts, totalReadCounts, wbcStatus,
    nSamplingEvents
  )


  dist_histogram <- parallel::mclapply(distance_statistics,
    FUN = function(input_list_elements) {
      return(input_list_elements[1])
    }
  ) %>%
    unlist()

  totalNumberOfSplits <- parallel::mclapply(distance_statistics,
    FUN = function(input_list_elements) {
      return(input_list_elements[2])
    }
  ) %>%
    unlist() %>%
    sum()

  StatisticsOfMutationPlacement <- parallel::mclapply(distance_statistics,
    FUN =
      function(input_list_elements) {
        return(input_list_elements[3])
      }
  ) %>%
    unlist()


  totalNumberOfSamplingEvents <- nSamplingEvents * length(postSampling)


  data <-
    data.frame(StatisticsOfMutationPlacement = StatisticsOfMutationPlacement)


  sum(is.na(data$StatisticsOfMutationPlacement))
  class(data$StatisticsOfMutationPlacement)


  ggplot(data = data, aes(x = StatisticsOfMutationPlacement, y = 1)) +
    geom_point()


  tryCatch(
    expr = {
      histo <-
        ggplot(data, aes(x = StatisticsOfMutationPlacement)) +
        geom_histogram(
          bins = 10, fill = "skyblue", color = "skyblue", alpha = 0.7
        ) +
        xlab("Splitting score") +
        ylab("total count") +
        ggtitle("Posterior sampling of branching probabilites") +
        geom_vline(
          xintercept = mean(StatisticsOfMutationPlacement),
          color = "blue", linetype = "dashed", linewidth = 1
        ) +
        labs(subtitle = sprintf("Tree %s - %s", treeName, clusterName)) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.subtitle = element_text(size = 18),
          axis.text = element_text(size = 16)
        )
      hist_data <- ggplot_build(histo)$data[[1]]
      max_y <- max(hist_data$count)
      histo <- histo + annotate("text",
        x = mean(StatisticsOfMutationPlacement) + 0.08,
        y = 0.9 * max_y, label = "mean", color = "blue",
        size = 7
      )
      print(histo)
    },
    error = function(e) {
      histo <- ggplot(data, aes(x = log(StatisticsOfMutationPlacement))) +
        geom_histogram(
          bins = 10, fill = "skyblue", color = "skyblue",
          alpha = 0.7
        ) +
        xlab("log(Splitting Score") +
        ylab("total count") +
        ggtitle("Posterior sampling of branching probabilites - Logarithmic
                Scale") +
        geom_vline(
          xintercept = log(mean(StatisticsOfMutationPlacement)),
          color = "blue", linetype = "dashed", linewidth = 1
        ) +
        labs(
          subtitle = sprintf("Tree %s - %s", treeName, clusterName),
          caption = "mean indicated by dashed red line"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.subtitle = element_text(size = 18),
          axis.text = element_text(size = 16)
        )
      hist_data <- ggplot_build(histo)$data[[1]]
      max_y <- max(hist_data$count)
      histo <- histo + annotate("text",
        x = log(mean(StatisticsOfMutationPlacement)) +
          0.08, y = 0.9 * max_y, label = "log(mean)",
        color = "blue", size = 7
      )
      print(histo)
    }
  )




  return(list(
    splittingFraction =
      totalNumberOfSplits / totalNumberOfSamplingEvents,
    branchingStatistics = StatisticsOfMutationPlacement
  ))
}




#' This function identifies cells that belong to the same CTC cluster - also
#' those which have been physically split. For each pair of tumour cells from
#' the same CTC cluster, the distnace postior is computed.
#'
#' @param sampleDescription A data frame with the description of each sample.
#' Expects the following columns:
#' Cluster: numeric vector indicating the cluster identity. Physically separated
#' clusters usually have different cluster identities, but this is not
#' necessary.
#' @param postSampling The loaded posterior sampling table.
#' @param treeName A string with the name of the tree that is output to the
#' plots.
#' @param nCells The total number of cells in the experiment.
#' @param nMutations The total number of mutations in the experiment.
#' @param nClusters The total number of clusters in the experiment.
#' @param alleleCount A numeric vector which indicates the number of alleles in
#' each of the clusters.
#' @param mutatedReadCounts A tibble containing the mutated reads. Rows are
#' mutations and columns are samples (clusters).
#' @param totalReadCounts A tibble containing the total read counts.
#' @param nMutationSamplingEvents The number of mutation that should be sampled
#' per tree.
#' @param nTreeSamplingEvents The number of trees that should be sampled.
#' @param cellPairSelection An optional parameter that takes a list of
#' pairs of strings-valued names of cells that should be analysed (the names as
#' in the samples_nodeDescription.tsv file).
#' It can also take a character vector, in which case the entries should be the
#' color coded names of the clusters.
#'
#' @return splittinProbs a vector that gives for each pair of cells the fraction
#' of trees for which they split
#' aggregatedBranchingProbabilities: a vector of aggregated probabilities for
#' all considered pairs of leaves and all sampled trees. At the moment only
#' implement if  cellPairSelection
#' parameter is passed to the function.
#' @export
#'
#' @examples
computeClusterSplits <- function(sampleDescription, postSampling, treeName,
                                 nCells, nMutations, nClusters, alleleCount,
                                 mutatedReadCounts, totalReadCounts,
                                 nMutationSamplingEvents = 1000,
                                 nTreeSamplingEvents = 500,
                                 cellPairSelection = NA) {
  desired_values <- sample(1:length(postSampling),
    size = nTreeSamplingEvents,
    replace = FALSE
  ) %>% sort()

  postSampling <- postSampling[desired_values]
  splittingProbs <- matrix(0, nrow = 0, ncol = 2) %>% as.data.frame()
  colnames(splittingProbs) <- c("Cluster", "Splitting_probability")
  aggregatedProbabilities <- vector()
  if (class(cellPairSelection) == "list") {
    counter <- 1
    system.time(
      for (it in cellPairSelection) {
        leaf1 <- which(sampleDescription$ClusterName == it[1]) - 1
        leaf2 <- which(sampleDescription$ClusterName == it[2]) - 1

        print(paste(paste("Computing genomic distances of leaves:", leaf1,
          sep = " "
        ), leaf2, sep = " "))
        posterior <- produce_Distance_Posterior(leaf1, leaf2, postSampling,
          treeName, nCells, nMutations,
          nClusters, alleleCount,
          sampleDescription$Cluster,
          mutatedReadCounts,
          totalReadCounts,
          sampleDescription$WBC,
          nSamplingEvents =
            nMutationSamplingEvents
        )
        splittingProbs <- rbind(
          splittingProbs,
          data.frame(
            Cluster = as.character(counter),
            Splitting_probability =
              posterior$splittingFraction
          )
        )
        aggregatedProbabilities <- c(
          aggregatedProbabilities,
          posterior$branchingStatistics
        )
        counter <- counter + 1
      }
    )
  } else if (class(cellPairSelection) == "character") {
    CTCclusters <- unique(cellPairSelection)
    CTCclusters <- CTCclusters[!(CTCclusters %in% c("ghostwhite", "gray93"))]
    system.time(
      for (it in CTCclusters) {
        cellsInCluster <- which(sampleDescription$color %in% it) - 1
        ## Make sure array indication is compatible with cpp
        cluster_done <- 0
        for (i in cellsInCluster) {
          if (cluster_done == 1) {
            cluster_done <- 0
            break
          }
          if (sampleDescription$WBC[i + 1] == 1) next
          j <- cellsInCluster[1]
          while (j < i) {
            if (cluster_done == 1) {
              break
            }
            if (sampleDescription$WBC[j + 1] == 1) {
              j <- j + 1
              next
            }
            print(paste(paste("Computing genomic distances of leaves:", i,
              sep = " "
            ), j, sep = " "))
            posterior <- produce_Distance_Posterior(i, j, postSampling,
              treeName, nCells,
              nMutations, nClusters,
              alleleCount,
              sampleDescription$Cluster,
              mutatedReadCounts,
              totalReadCounts,
              sampleDescription$WBC,
              nSamplingEvents =
                nMutationSamplingEvents,
              clusterName = it
            )
            splittingProbs <- rbind(
              splittingProbs,
              data.frame(
                Cluster = it,
                Splitting_probability =
                  posterior$splittingFraction
              )
            )
            aggregatedProbabilities <- c(
              aggregatedProbabilities,
              posterior$branchingStatistics
            )
            j <- j + 1
            cluster_done <- 1
          }
        }
      }
    )
  } else {
    CTCclusters <- unique(sampleDescription$color)
    CTCclusters <- CTCclusters[!(CTCclusters %in% c("ghostwhite", "gray93"))]
    system.time(
      for (it in CTCclusters) {
        cellsInCluster <- which(sampleDescription$color %in% it) - 1
        ## Make sure array indication is compatible with cpp
        # cluster_done <- 0
        for (i in cellsInCluster) {
          # if(cluster_done == 1){
          #  cluster_done <- 0
          #  break
          # }
          if (sampleDescription$WBC[i + 1] == 1) next
          j <- cellsInCluster[1]
          while (j < i) {
            # if(cluster_done == 1){
            #  break
            # }
            if (sampleDescription$WBC[j + 1] == 1) {
              j <- j + 1
              next
            }
            print(paste(paste("Computing genomic distances of leaves:", i,
              sep = " "
            ), j, sep = " "))
            posterior <- produce_Distance_Posterior(i, j, postSampling,
              treeName, nCells,
              nMutations, nClusters,
              alleleCount,
              sampleDescription$Cluster,
              mutatedReadCounts,
              totalReadCounts,
              sampleDescription$WBC,
              nSamplingEvents =
                nMutationSamplingEvents,
              clusterName = it
            )
            splittingProbs <- rbind(
              splittingProbs,
              data.frame(
                Cluster = it,
                Splitting_probability =
                  posterior$splittingFraction
              )
            )
            j <- j + 1
            # cluster_done <- 1
          }
        }
      }
    )
  }



  return(list(
    splittingProbs = splittingProbs,
    aggregatedBranchingProbabilities = aggregatedProbabilities
  ))
}




#' Loads all necessary data for the CTC-project.
#' Specifically it return a named list as follows:
#' postSampling: Loads the posterior sampling tsv as a list of named vectors
#' with the following columns: the (unnormalised) LogScore, estimated sequencing
#' error rate, the estimated dropout rate, logTau and the Tree in parent vector
#' format meaning that the i'th entry of the vector is te parent node of the
#' entry i.
#' Nodes are counted from zero and the root is length(Tree)
#'
#' @param inputFolder The total number of CTC-clusters
#' @param treeName
#'
#' @return postSampling: Loads the posterior sampling tsv as a list of named
#' vectors with the following columns: the (unnormalised) LogScore, estimated
#' sequencing error rate, the estimated dropout rate, logTau and the Tree in
#' parent vector format meaning that the i'th entry of the vector is the parent
#' node of the entry i. Nodes are counted from zero and the root is length(Tree)
#' @export
#'
#' @examples
load_data <- function(inputFolder, treeName) {
  ## Define paths

  posteriorSamplingFile <- sprintf(
    "%s/%s/%s_postSampling.tsv", inputFolder,
    treeName, treeName
  )

  countFile <- sprintf("%s/%s/%s.txt", inputFolder, treeName, treeName)
  descriptionFile <- sprintf(
    "%s/%s/%s_samples_nodeDescription.tsv",
    inputFolder, treeName, treeName
  )


  ## Load data

  postSampling <- read_delim(posteriorSamplingFile,
    delim = "\t", col_names = c(
      "LogScore", "SequencingErrorRate",
      "DropoutRate", "LogTau", "Tree"
    )
  )
  postSampling <- split(postSampling, seq(nrow(postSampling)))


  counts <- read_delim(countFile,
    delim = "\t", col_names = FALSE
  )
  description <- read_delim(descriptionFile,
    delim = "\t", col_names = c(
      "Cluster", "CellCount", "TCs", "WBCs",
      "Description"
    )
  )
  nCells <- sum(description$CellCount)
  nClusters <- nrow(description)
  nMutations <- nrow(counts)
  alleleCount <- description$CellCount * 2


  description <- description %>%
    mutate(color = regmatches(Description, regexpr(
      "color=([a-zA-Z]+[0-9]*)",
      Description
    )) %>%
      substr(start = 7, stop = (nchar(.))))



  ClusterID <- vector()
  for (i in 1:nClusters) {
    ClusterID <- c(
      ClusterID,
      rep.int(
        i - 1,
        description$CellCount[i]
      )
    )
  }
  ## Note that Cpp counts arrays from zero, so the cluster IDs are counted
  ## likewise in order to be compatible with Cpp code.

  ## Pull apart the count file into counts for mutated read and total counts
  ## respectively
  mutatedReadCounts <- matrix(0, nrow = nMutations, ncol = 0)
  for (j in 1:nClusters) {
    mutatedReadCounts <- cbind(mutatedReadCounts, counts[, 4 + 2 * j])
  }

  totalReadCounts <- matrix(0, nrow = nMutations, ncol = 0)
  for (j in 1:nClusters) {
    totalReadCounts <- cbind(totalReadCounts, counts[, 4 + 2 * j - 1])
  }


  wildtypeReadCounts <- totalReadCounts - mutatedReadCounts


  mutatedReadCounts <- mutatedReadCounts %>%
    t() %>%
    as.data.frame() %>%
    as.list()
  wildtypeReadCounts <- wildtypeReadCounts %>%
    t() %>%
    as.data.frame() %>%
    as.list()
  totalReadCounts <- totalReadCounts %>%
    t() %>%
    as.data.frame() %>%
    as.list()


  mutationDescription <- counts[, 1:4]

  ## wbc status indicates which of the cells is a white blood cells and which
  ## one isn't.
  ## So far, the cells are arbitrary, and I will assign the fist cells from a
  ## cluster to be WBCs.
  wbcStatus <- rep(0, nCells)

  for (i in 1:nClusters) {
    j <- 1
    while (j <= description$WBCs[i]) { # Iterating over the number of White
      # blood cells of a cluster
      wbcStatus[which(ClusterID == i - 1)[1] + j - 1] <- 1 # and identifying the
      # first cell
      # that belongs to a cluster and counting from then on
      ## Note: The cluster IDs are counted from zero!
      j <- j + 1
    }
  }



  sample_description <- data.frame(
    Cluster = ClusterID,
    ClusterName = description$Cluster[ClusterID + 1],
    WBC = wbcStatus,
    color = description$color[ClusterID + 1]
  )

  sample_description <- sample_description %>%
    mutate(single_cell = !(duplicated(Cluster)) &
      !(duplicated(Cluster, fromLast = TRUE)))




  return(list(
    "postSampling" = postSampling, "nClusters" = nClusters,
    "clusterID" = ClusterID, "nCells" = nCells,
    "nMutations" = nMutations, "nClusters" = nClusters,
    "alleleCount" = alleleCount,
    "mutatedReadCounts" = mutatedReadCounts,
    "totalReadCounts" = totalReadCounts, "wbcStatus" = wbcStatus,
    "sample_description" = sample_description,
    "mutationDescription" = mutationDescription,
    #              "annotations" = annotations,
    "sampleName" = treeName, "directory" = inputFolder
  ))
}







#' Takes called genotypes in .ped format, computes a pairwise distance matrix
#' and indentifies pairs of distinct cells (or cell clusters, needs a manual
#' check) that are genetically similar to each other. Similar means that their
#' genetic distance lies in the 1% quantile of the set of all pairwise genetic
#' distances.
#' As the distance the Hamming distance is chosen.
#'
#' @param inputFolder
#' @param treeName
#'
#' @return
#' monoclonal_pairs: A list of pairs of cell names that are similar to each
#' other.
#' distance_matrix: A matrix indicates all pairwise distnaces of suggested
#' genotypes.
#' full_distance_matrix: The full pairwise distance matrix of all genotypes.
#'
#'
#' @export
#'
#' @examples
load_monoclonal_pairs <- function(inputFolder, treeName, cutoff = "") {
  data_file <- sprintf(
    "%s/%s/%s_genotypes.ped", inputFolder, treeName,
    treeName
  )

  data <- read_delim(data_file, delim = "\t", col_names = FALSE)

  data2 <- data %>% select(!2:6)

  distance_matrix <- matrix(0, nrow = nrow(data2), ncol = nrow(data2))


  for (i in 1:nrow(data2)) {
    j <- 1
    while (j < i) {
      row_i <- data2 %>%
        select(!1) %>%
        slice(i)
      row_j <- data2 %>%
        select(!1) %>%
        slice(j)

      distance_matrix[i, j] <- sum(!(row_i == row_j))
      j <- j + 1
    }
  }


  distance_vector <- as.vector(distance_matrix[lower.tri(distance_matrix)])



  if (class(cutoff) == "numeric") {
    monoclonal_candidate_cutoff <- cutoff
  } else {
    monoclonal_candidate_cutoff <- quantile(distance_vector, probs = 0.01)
  }


  sum(distance_vector <= monoclonal_candidate_cutoff)
  which(distance_vector <= monoclonal_candidate_cutoff)

  print("1% quantile of genetic distances:")
  print(monoclonal_candidate_cutoff)

  plot(
    ggplot(data.frame(x = distance_vector), aes(x = x)) +
      geom_histogram(binwidth = 2) +
      geom_vline(
        xintercept = monoclonal_candidate_cutoff, linetype = "dashed",
        color = "red"
      )
  )


  candidates <- list()
  candidate_index <- vector()
  iterator <- 0
  number_of_output_pairs <- 15
  for (count in 0:monoclonal_candidate_cutoff) {
    all_elements <- which(distance_matrix == count)
    all_elements_list <- list()
    for (it in all_elements) {
      coordinates1 <- ((it - 1) %% dim(distance_matrix)[2]) + 1
      coordinates2 <- ((it - 1) %/% dim(distance_matrix)[2]) + 1
      all_elements_list <-
        append(all_elements_list, list(c(coordinates1, coordinates2)))
    }

    for (it in all_elements_list) {
      if (it[1] <= it[2]) next


      # Check whether the candidate pair of cells consists of single tumour
      # cells:



      candidates <-
        c(candidates, list(c(
          as.character(data2[it[1], 1]),
          as.character(data2[it[2], 1])
        )))
      candidate_index <- c(candidate_index, it[1], it[2])

      iterator <- iterator + 1

      if (iterator == number_of_output_pairs) break
    }
    if (iterator == number_of_output_pairs) break
  }
  if (length(unique(sort(candidate_index))) != 0) {
    distance_matrix2 <-
      distance_matrix[
        unique(sort(candidate_index)),
        unique(sort(candidate_index))
      ]
    colnames(distance_matrix2) <- data2[unique(sort(candidate_index)), 1]$X1
  } else {
    distance_matrix2 <- 0
  }


  distance_matrix <- as.data.frame(distance_matrix)
  colnames(distance_matrix) <- data2$X1
  rownames(distance_matrix) <- data2$X1
  return(
    list(
      monoclonal_pairs = candidates,
      distance_matrix = distance_matrix2,
      full_distance_matrix = distance_matrix
    )
  )
}
