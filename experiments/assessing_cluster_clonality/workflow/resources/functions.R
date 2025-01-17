############
# Function Definitions
############

library(Rcpp)
library(tidyverse)

sourceCpp("~/work/CTC-SCITE/experiments/assessing_cluster_clonality/workflow/resources/mutations_placement.cpp")



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
is_driver <- function(mutations, annotations) {
  annotated_mutations <- annotations |>
    filter(
      annotations$"#CHROM" == as.character(mutations[1]) &
        annotations$POS == as.numeric(mutations[2])
    )

  check <- annotated_mutations |>
    dplyr::select(c(
      "CGI-Oncogenic Summary", "CGI-Oncogenic Prediction",
      "CGI-External oncogenic annotation"
    )) %in%
    c("oncogenic (predicted)", "driver (oncodriveMUT)") |>
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
#' @param post_sampling loaded list of tibbles tibble containing the posterior
#'                     Sampling
#' @param tree_name character string: Name of the tree for the output plot
#' @param n_cells integer-valued total number of cells in the dataset
#' @param n_mutations integer-valued total number of mutations in the dataset
#' @param n_clusters integer-valued total number of clusters in the dataset
#' @param allele_count integer vector of numbers of alleles per clusters
#' @param cluster_id integer vector of cluster IDs
#' @param mutated_read_counts list of integer-valued vectors indicating the
#' number of mutated read per mutation (list index) and sample (vector index)
#' @param total_read_counts list of integer-valued vectors indicating the total
#' number of reads per mutation (list index) and sample (vector index)
#' @param wbc_status boolean vector of length n_cells indicating for each cell
#' if it is a white blood cell (TRUE) or not (FALSE)
#'
#' @return splittingFraction: The fraction of sampling events for which the pair
#'         of cells
#' shows branching evolution
#'
#' @export
#'
#' @examples
produce_distance_posterior <- function(leaf1, leaf2, post_sampling, tree_name,
                                       n_cells, n_mutations, n_clusters,
                                       allele_count, cluster_id,
                                       mutated_read_counts, total_read_counts,
                                       wbc_status, n_sampling_events = 20,
                                       cluster_name = "") {
  ## For each row in the posterior Sampling file, the distance of two leaves is
  ## computed

  print("Computing the posterior distribution")

  distance_statistics <- parallel::mclapply(post_sampling,
    FUN = computePairwiseDistanceOfLeavesGivenTree, leaf1, leaf2,
    n_cells, n_mutations, n_clusters, allele_count,
    cluster_id, mutated_read_counts, total_read_counts, wbc_status,
    n_sampling_events
  )


  no_splits <- parallel::mclapply(distance_statistics,
    FUN = function(input_list_elements) {
      return(input_list_elements[2])
    }
  ) |>
    unlist() |>
    sum()

  mutation_placement_stats <- parallel::mclapply(distance_statistics,
    FUN =
      function(input_list_elements) {
        return(input_list_elements[3])
      }
  ) |>
    unlist()


  no_sampling_events <- n_sampling_events * length(post_sampling)


  data <-
    data.frame(
      mutation_placement_stats = mutation_placement_stats
    )


  sum(is.na(data$mutation_placement_stats))
  class(data$mutation_placement_stats)


  ggplot2::ggplot(
    data = data, ggplot2::aes(x = mutation_placement_stats, y = 1)
  ) +
    ggplot2::geom_point()


  tryCatch(
    expr = {
      histo <-
        ggplot2::ggplot(
          data, ggplot2::aes(x = mutation_placement_stats)
        ) +
        ggplot2::geom_histogram(
          bins = 10, fill = "skyblue", color = "skyblue", alpha = 0.7
        ) +
        ggplot2::xlab("Splitting score") +
        ggplot2::ylab("total count") +
        ggplot2::ggtitle("Posterior sampling of branching probabilites") +
        ggplot2::geom_vline(
          xintercept = mean(mutation_placement_stats),
          color = "blue", linetype = "dashed", linewidth = 1
        ) +
        ggplot2::labs(
          subtitle = sprintf("Tree %s - %s", tree_name, cluster_name)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20, face = "bold"),
          axis.title.x = ggplot2::element_text(size = 18),
          axis.title.y = ggplot2::element_text(size = 18),
          plot.subtitle = ggplot2::element_text(size = 18),
          axis.text = ggplot2::element_text(size = 16)
        )
      hist_data <- ggplot2::ggplot_build(histo)$data[[1]]
      max_y <- max(hist_data$count)
      histo <- histo + ggplot2::annotate("text",
        x = mean(mutation_placement_stats) + 0.08,
        y = 0.9 * max_y, label = "mean", color = "blue",
        size = 7
      )
      print(histo)
    },
    error = function(e) {
      histo <- ggplot2::ggplot(
        data, ggplot2::aes(x = log(mutation_placement_stats))
      ) +
        ggplot2::geom_histogram(
          bins = 10, fill = "skyblue", color = "skyblue",
          alpha = 0.7
        ) +
        ggplot2::xlab("log(Splitting Score") +
        ggplot2::ylab("total count") +
        ggplot2::ggtitle(
          "Posterior sampling of branching probabilites - Logarithmic Scale"
        ) +
        ggplot2::geom_vline(
          xintercept = log(mean(mutation_placement_stats)),
          color = "blue", linetype = "dashed", linewidth = 1
        ) +
        ggplot2::labs(
          subtitle = sprintf("Tree %s - %s", tree_name, cluster_name),
          caption = "mean indicated by dashed red line"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 20, face = "bold"),
          axis.title.x = ggplot2::element_text(size = 18),
          axis.title.y = ggplot2::element_text(size = 18),
          plot.subtitle = ggplot2::element_text(size = 18),
          axis.text = ggplot2::element_text(size = 16)
        )
      hist_data <- ggplot2::ggplot_build(histo)$data[[1]]
      max_y <- max(hist_data$count)
      histo <- histo + ggplot2::annotate("text",
        x = log(mean(mutation_placement_stats)) +
          0.08, y = 0.9 * max_y, label = "log(mean)",
        color = "blue", size = 7
      )
      print(histo)
    }
  )




  return(list(
    splittingFraction =
      no_splits / no_sampling_events,
    branchingStatistics = mutation_placement_stats
  ))
}




#' This function identifies cells that belong to the same CTC cluster - also
#' those which have been physically split. For each pair of tumour cells from
#' the same CTC cluster, the distnace postior is computed.
#'
#' @param sample_description A data frame with the description of each sample.
#' Expects the following columns:
#' Cluster: numeric vector indicating the cluster identity. Physically separated
#' clusters usually have different cluster identities, but this is not
#' necessary.
#' @param post_sampling The loaded posterior sampling table.
#' @param tree_name A string with the name of the tree that is output to the
#' plots.
#' @param n_cells The total number of cells in the experiment.
#' @param n_mutations The total number of mutations in the experiment.
#' @param n_clusters The total number of clusters in the experiment.
#' @param allele_count A numeric vector which indicates the number of alleles in
#' each of the clusters.
#' @param mutated_read_counts A tibble containing the mutated reads. Rows are
#' mutations and columns are samples (clusters).
#' @param total_read_counts A tibble containing the total read counts.
#' @param n_mutation_sampling_events The number of mutation that should be
#' sampled per tree.
#' @param n_tree_sampling_events The number of trees that should be sampled.
#' @param cell_pair_selection An optional parameter that takes a list of
#' pairs of strings-valued names of cells that should be analysed (the names as
#' in the samples_nodeDescription.tsv file).
#' It can also take a character vector, in which case the entries should be the
#' color coded names of the clusters.
#'
#' @return splittinProbs a vector that gives for each pair of cells the fraction
#' of trees for which they split
#' aggregatedBranchingProbabilities: a vector of aggregated probabilities for
#' all considered pairs of leaves and all sampled trees. At the moment only
#' implement if  cell_pair_selection
#' parameter is passed to the function.
#' @export
#'
#' @examples
compute_cluster_splits <- function(sample_description, post_sampling, tree_name,
                                   n_cells, n_mutations, n_clusters,
                                   allele_count, mutated_read_counts,
                                   total_read_counts,
                                   n_mutation_sampling_events = 1000,
                                   n_tree_sampling_events = 500,
                                   cell_pair_selection = NA) {
  desired_values <- sample(seq_along(post_sampling),
    size = n_tree_sampling_events,
    replace = FALSE
  ) |> sort()

  post_sampling <- post_sampling[desired_values]
  splitting_probs <- matrix(0, nrow = 0, ncol = 2) |> as.data.frame()
  colnames(splitting_probs) <- c("Cluster", "Splitting_probability")
  aggregated_probabilities <- vector()
  if (class(cell_pair_selection) == "list") {
    counter <- 1
    system.time(
      for (it in cell_pair_selection) {
        print(it)
        leaf1 <- which(sample_description$ClusterName == it[1]) - 1
        leaf2 <- which(sample_description$ClusterName == it[2]) - 1

        print(
          paste(
            "Computing genomic distances of leaves:", leaf1, leaf2,
            sep = " "
          )
        )
        posterior <- produce_distance_posterior(leaf1, leaf2, post_sampling,
          tree_name, n_cells, n_mutations,
          n_clusters, allele_count,
          sample_description$Cluster,
          mutated_read_counts,
          total_read_counts,
          sample_description$WBC,
          n_sampling_events =
            n_mutation_sampling_events
        )
        splitting_probs <- rbind(
          splitting_probs,
          data.frame(
            Cluster = as.character(counter),
            Splitting_probability =
              posterior$splittingFraction
          )
        )
        aggregated_probabilities <- c(
          aggregated_probabilities,
          posterior$branchingStatistics
        )
        counter <- counter + 1
      }
    )
  } else if (class(cell_pair_selection) == "character") {
    ctc_clusters <- unique(cell_pair_selection)
    ctc_clusters <- ctc_clusters[!(ctc_clusters %in% c("ghostwhite", "gray93"))]
    print(ctc_clusters)

    system.time(
      for (it in ctc_clusters) {
        cells_in_cluster <- which(sample_description$color == it) - 1
        ## Make sure array indication is compatible with cpp
        cluster_done <- 0
        for (i in cells_in_cluster) {
          if (cluster_done == 1) {
            cluster_done <- 0
            break
          }
          if (sample_description$WBC[i + 1] == 1) next
          j <- cells_in_cluster[1]
          while (j < i) {
            if (cluster_done == 1) {
              break
            }
            if (sample_description$WBC[j + 1] == 1) {
              j <- j + 1
              next
            }
            print(
              paste("Computing genomic distances of leaves:", i, j, sep = " ")
            )

            posterior <- produce_distance_posterior(i, j, post_sampling,
              tree_name, n_cells,
              n_mutations, n_clusters,
              allele_count,
              sample_description$Cluster,
              mutated_read_counts,
              total_read_counts,
              sample_description$WBC,
              n_sampling_events =
                n_mutation_sampling_events,
              cluster_name = it
            )

            splitting_probs <- rbind(
              splitting_probs,
              data.frame(
                Cluster = it,
                Splitting_probability =
                  posterior$splittingFraction
              )
            )
            aggregated_probabilities <- c(
              aggregated_probabilities,
              posterior$branchingStatistics
            )
            j <- j + 1
            cluster_done <- 1
          }
        }
      }
    )
  } else {
    ctc_clusters <- unique(sample_description$color)
    ctc_clusters <- ctc_clusters[!(ctc_clusters %in% c("ghostwhite", "gray93"))]
    system.time(
      for (it in ctc_clusters) {
        cells_in_cluster <- which(sample_description$color %in% it) - 1
        ## Make sure array indication is compatible with cpp

        for (i in cells_in_cluster) {
          if (sample_description$WBC[i + 1] == 1) next
          j <- cells_in_cluster[1]
          while (j < i) {
            if (sample_description$WBC[j + 1] == 1) {
              j <- j + 1
              next
            }
            print(
              paste("Computing genomic distances of leaves:", i, j, sep = " ")
            )
            posterior <- produce_distance_posterior(i, j, post_sampling,
              tree_name, n_cells,
              n_mutations, n_clusters,
              allele_count,
              sample_description$Cluster,
              mutated_read_counts,
              total_read_counts,
              sample_description$WBC,
              n_sampling_events =
                n_mutation_sampling_events,
              cluster_name = it
            )
            print("Posterior computed")
            splitting_probs <- rbind(
              splitting_probs,
              data.frame(
                Cluster = it,
                Splitting_probability =
                  posterior$splittingFraction
              )
            )
            j <- j + 1
          }
        }
      }
    )
  }



  return(list(
    splitting_probs = splitting_probs,
    aggregatedBranchingProbabilities = aggregated_probabilities
  ))
}




#' Loads all necessary data for the CTC-project.
#' Specifically it return a named list as follows:
#' post_sampling: Loads the posterior sampling tsv as a list of named vectors
#' with the following columns: the (unnormalised) LogScore, estimated sequencing
#' error rate, the estimated dropout rate, logTau and the Tree in parent vector
#' format meaning that the i'th entry of the vector is te parent node of the
#' entry i.
#' Nodes are counted from zero and the root is length(Tree)
#'
#' @param input_folder The total number of CTC-clusters
#' @param tree_name
#'
#' @return post_sampling: Loads the posterior sampling tsv as a list of named
#' vectors with the following columns: the (unnormalised) LogScore, estimated
#' sequencing error rate, the estimated dropout rate, logTau and the Tree in
#' parent vector format meaning that the i'th entry of the vector is the parent
#' node of the entry i. Nodes are counted from zero and the root is length(Tree)
#' @export
#'
#' @examples
load_data <- function(input_folder, tree_name) {
  ## Define paths

  posterior_sampling_file <- sprintf(
    "%s/%s/%s_postSampling.tsv", input_folder,
    tree_name, tree_name
  )

  count_file <- sprintf("%s/%s/%s.txt", input_folder, tree_name, tree_name)
  description_file <- sprintf(
    "%s/%s/%s_samples_nodeDescription.tsv",
    input_folder, tree_name, tree_name
  )


  ## Load data

  post_sampling <- readr::read_delim(posterior_sampling_file,
    delim = "\t", col_names = c(
      "LogScore", "SequencingErrorRate",
      "DropoutRate", "LogTau", "Tree"
    )
  )
  post_sampling <- split(post_sampling, seq_len(nrow(post_sampling)))


  counts <- readr::read_delim(count_file,
    delim = "\t", col_names = FALSE
  )
  description <- readr::read_delim(description_file,
    delim = "\t", col_names = c(
      "Cluster", "CellCount", "TCs", "WBCs",
      "Description"
    )
  )

  n_cells <- sum(description$CellCount)
  n_clusters <- nrow(description)
  n_mutations <- nrow(counts)
  allele_count <- description$CellCount * 2


  description <- description %>%
    dplyr::mutate(color = regmatches(Description, regexpr(
      "color=([a-zA-Z]+[0-9]*)",
      Description
    )) %>%
    substr(start = 7, stop = (nchar(.))))


  cluster_id <- vector()
  for (i in seq_len(n_clusters)) {
    cluster_id <- c(
      cluster_id,
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
  mutated_read_counts <- matrix(0, nrow = n_mutations, ncol = 0)
  for (j in seq_len(n_clusters)) {
    mutated_read_counts <- cbind(mutated_read_counts, counts[, 4 + 2 * j])
  }

  total_read_counts <- matrix(0, nrow = n_mutations, ncol = 0)
  for (j in seq_len(n_clusters)) {
    total_read_counts <- cbind(total_read_counts, counts[, 4 + 2 * j - 1])
  }


  wildtype_read_counts <- total_read_counts - mutated_read_counts


  mutated_read_counts <- mutated_read_counts |>
    t() |>
    as.data.frame() |>
    as.list()
  wildtype_read_counts <- wildtype_read_counts |>
    t() |>
    as.data.frame() |>
    as.list()
  total_read_counts <- total_read_counts |>
    t() |>
    as.data.frame() |>
    as.list()


  mutation_description <- counts[, 1:4]

  ## wbc status indicates which of the cells is a white blood cells and which
  ## one isn't.
  ## So far, the cells are arbitrary, and I will assign the fist cells from a
  ## cluster to be WBCs.
  wbc_status <- rep(0, n_cells)

  for (i in seq_len(n_clusters)) {
    j <- 1
    while (j <= description$WBCs[i]) { # Iterating over the number of White
      # blood cells of a cluster
      wbc_status[which(cluster_id == i - 1)[1] + j - 1] <- 1
      # and identifying the first cell that belongs to a cluster and counting
      # from then on.
      ## Note: The cluster IDs are counted from zero!
      j <- j + 1
    }
  }



  sample_description <- data.frame(
    Cluster = cluster_id,
    ClusterName = description$Cluster[cluster_id + 1],
    WBC = wbc_status,
    color = description$color[cluster_id + 1]
  )

  sample_description <- sample_description |>
    dplyr::mutate(
      single_cell =
        !(duplicated(Cluster)) &
          !(duplicated(Cluster, fromLast = TRUE))
    )




  return(list(
    "post_sampling" = post_sampling, "n_clusters" = n_clusters,
    "cluster_id" = cluster_id, "n_cells" = n_cells,
    "n_mutations" = n_mutations, "allele_count" = allele_count,
    "mutated_read_counts" = mutated_read_counts,
    "total_read_counts" = total_read_counts, "wbc_status" = wbc_status,
    "sample_description" = sample_description,
    "mutation_description" = mutation_description,
    "sampleName" = tree_name, "directory" = input_folder
  ))
}







#' Takes called genotypes in .ped format, computes a pairwise distance matrix
#' and indentifies pairs of distinct cells (or cell clusters, needs a manual
#' check) that are genetically similar to each other. Similar means that their
#' genetic distance lies in the 1% quantile of the set of all pairwise genetic
#' distances.
#' As the distance the Hamming distance is chosen.
#'
#' @param input_folder
#' @param tree_name
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
load_monoclonal_pairs <- function(input_folder, tree_name, cutoff = "") {
  data_file <- sprintf(
    "%s/%s/%s_genotypes.ped", input_folder, tree_name,
    tree_name
  )

  data <- readr::read_delim(data_file, delim = "\t", col_names = FALSE)

  data2 <- data |> dplyr::select(!2:6)

  distance_matrix <- matrix(0, nrow = nrow(data2), ncol = nrow(data2))


  for (i in seq_len(nrow(data2))) {
    j <- 1
    while (j < i) {
      row_i <- data2 |>
        dplyr::select(!1) |>
        dplyr::slice(i)
      row_j <- data2 |>
        dplyr::select(!1) |>
        dplyr::slice(j)

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
    ggplot2::ggplot(
      data.frame(x = distance_vector), ggplot2::aes(x = rlang::.data$x)
    ) +
      ggplot2::geom_histogram(binwidth = 2) +
      ggplot2::geom_vline(
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
