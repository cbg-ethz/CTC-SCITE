input_folder <- "/Users/jgawron/Documents/projects/CTC_backup/input_folder"

simulation_input_folder <- "/Users/jgawron/Documents/projects/CTC_backup/simulations/simulations2"
tree_name <- "Br16_AC"
n_sampling_events <- 100

source("functions.R")


tree_names <- c("Br11", "Br16_AC", "Br16_B", "Br16_C", "Br23", "Br26", "Br30", "Br37", "Br38", "Br39", "Br44", "Br45", "Br46", "Br53", "Br57", "Br61", "Br7", "Brx50", "LM2", "Lu2", "Lu7", "Ov8", "Pr6", "Pr9")


mean_splitting_scores_mono <- vector()
for (tree_name in tree_names) {
  # input <- load_data(input_folder, tree_name)

  # description_data <-
  #  read_delim(
  #    file.path(
  #      input$directory,
  ###      input$sampleName,
  #      paste0(input$sampleName, "_samples_nodeDescription.tsv")
  #    ),
  #    delim = "\t",
  ###    col_names = FALSE,
  #    quote = "none"
  #  )
  # colnames(description_data) <-
  #  c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description")


  folders <- list.dirs(path = simulation_input_folder, recursive = FALSE)
  matching_folders <- folders[grepl(tree_name, basename(folders))]

  for (idx1 in 1:length(matching_folders)) {
    simulation_instance <- basename(matching_folders[idx1])
    input_simulated <- load_data(simulation_input_folder, simulation_instance)
    sample_description_simulated <- input_simulated$sample_description
    for (color in c(
      "orchid", "orchid1", "orchid2",
      "orchid3", "orchid4", "darkorchid",
      "darkorchid1", "darkorchid2", "darkorchid3",
      "darkorchid4", "purple", "purple1",
      "purple2", "purple3", "purple4"
    )) {
      distance_simulated <-
        computeClusterSplits(
          input_simulated$sample_description, input_simulated$postSampling,
          simulation_instance, input_simulated$nCells, input_simulated$nMutations,
          input_simulated$nClusters, input_simulated$alleleCount,
          input_simulated$mutatedReadCounts, input_simulated$totalReadCounts,
          nMutationSamplingEvents = n_sampling_events,
          nTreeSamplingEvents = n_sampling_events,
          cellPairSelection = c(color)
        )
      plot(
        ggplot(
          data.frame(x = distance_simulated$aggregatedBranchingProbabilities), aes(x = x)
        ) +
          geom_histogram(binwidth = 0.01)
      )

      mean_splitting_scores_mono <- c(mean_splitting_scores_mono, mean(distance_simulated$aggregatedBranchingProbabilities))
    }
  }
}



involved_cell_indices <- sub(paste0(".*", tree_name, "_"), "", simulation_instance) %>%
  strsplit("_") %>%
  unlist()
involved_cell_indices <- as.numeric(involved_cell_indices)

involved_single_tumor_cells <- description_data$sample_name[involved_cell_indices]

pairs <- combn(involved_single_tumor_cells, 2, simplify = FALSE)

distance_separate <-
  computeClusterSplits(input$sample_description, input$postSampling,
    treeName, input$nCells, input$nMutations,
    input$nClusters, input$alleleCount,
    input$mutatedReadCounts, input$totalReadCounts,
    nMutationSamplingEvents = n_sampling_events,
    nTreeSamplingEvents = n_sampling_events,
    cellPairSelection = pairs
  )

plot(
  ggplot(
    data.frame(x = distance_separate$aggregatedBranchingProbabilities), aes(x = x)
  ) +
    geom_histogram(binwidth = 0.01)
)
mean(distance_separate$aggregatedBranchingProbabilities)




#########

mean_splitting_scores_mono <- mean_splitting_scores_mono[!is.na(mean_splitting_scores_mono)]

data.frame(splitting_probs = mean_splitting_scores_mono, Monoclonal = TRUE) %>%
  ggplot(aes(y = splitting_probs)) +
  geom_boxplot() +
  theme_minimal()

load("~/Documents/projects/CTC_backup/simulations/simulation3/mean_branching_probs_oligo.RData")
load("~/Documents/projects/CTC_backup/simulations/simulation3/mean_branching_probs_mono.RData")
load("~/Documents/projects/CTC_backup/simulations/simulation3/deviance_splitting_score.RData")

splitting_probs <- data.frame(splitting_probs = mean_splitting_scores_mono, Oligoclonal = "Monoclonal")
splitting_probs2 <- data.frame(splitting_probs = mean_branching_probs_oligo, Oligoclonal = "Oligoclonal")
splitting_probs_single_cells <- data.frame(splitting_probs = -deviance_splitting_score + mean_branching_probs_oligo, Oligoclonal = "Genetically distinct single cells")

splitting_probs <- rbind(splitting_probs, splitting_probs2, splitting_probs_single_cells)
splitting_probs <- splitting_probs %>% mutate(Oligoclonal = factor(Oligoclonal, levels = c("Monoclonal", "Oligoclonal", "Genetically distinct single cells")))

splitting_probs %>%
  ggplot(aes(y = splitting_probs, x = Oligoclonal, group = Oligoclonal)) +
  geom_boxplot() +
  ylab("Mean splitting probability") +
  xlab("Clonality status of CTC cluster") +
  theme_minimal()
theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
  axis.title.x = element_text(size = 20),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 18)
)

data.frame(y = deviance_splitting_score) %>%
  ggplot(aes(y = y)) +
  geom_boxplot() +
  theme_minimal()


data.frame(y = -deviance_splitting_score + mean_branching_probs_oligo) %>%
  ggplot(aes(y = y)) +
  geom_boxplot() +
  theme_minimal()
