input_folder <- "/Users/jgawron/Documents/projects/CTC_backup/input_folder"


tree_name <- "Br16_AC"
n_sampling_events <- 100

source("/Users/jgawron/Documents/projects/CTC-SCITE/experiments/assessing_cluster_clonality/workflow/resources/functions.R")


tree_names <- c("Pr9", "Br11", "Br16_AC", "Br16_B", "Br16_C", "Br23", "Br26", "Br38", "Br39", "Br57", "Br61", "Br7", "Brx50", "LM2", "Pr6")
folders <- list.dirs(path = simulation_input_folder, recursive = FALSE)

mean_splitting_scores_mono <- vector()
mean_splitting_scores_oligo <- vector()
deviance_splitting_score <- vector()

simulation_input_folder <- "/Users/jgawron/Documents/projects/CTC_backup/simulations/simulations2"
for (tree_name in tree_names) {
  print(tree_name)

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
        compute_cluster_splits(
          input_simulated$sample_description, input_simulated$post_sampling,
          simulation_instance, input_simulated$n_cells, input_simulated$n_mutations,
          input_simulated$n_clusters, input_simulated$allele_count,
          input_simulated$mutated_read_counts, input_simulated$total_read_counts,
          n_mutation_sampling_events = n_sampling_events,
          n_tree_sampling_events = n_sampling_events,
          cell_pair_selection = c(color)
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

mean_splitting_scores_mono <- mean_splitting_scores_mono[!is.na(mean_splitting_scores_mono)]

save(mean_splitting_scores_mono, file = "~/Documents/projects/CTC_backup/simulations/simulation3/mean_branching_probs_mono.RData")









simulation_input_folder <- "/Users/jgawron/Documents/projects/CTC_backup/simulations/simulations3_new"
for (tree_name in tree_names) {
  print(tree_name)

  input <- load_data(input_folder, tree_name)
  description_data <-
    read_delim(
      file.path(
        input$directory,
        input$sampleName,
        paste0(input$sampleName, "_samples_nodeDescription.tsv")
      ),
      delim = "\t",
      col_names = FALSE,
      quote = "none"
    )
  colnames(description_data) <-
    c("sample_name", "total_number_cells", "tumor_cells", "WBCs", "description")


  matching_folders <- folders[grepl(tree_name, basename(folders))]

  for (idx1 in 1:length(matching_folders)) {
    tryCatch(
      {
        simulation_instance <- basename(matching_folders[idx1])
        input_simulated <- load_data(simulation_input_folder, simulation_instance)
        sample_description_simulated <- input_simulated$sample_description
        # for (color in c(
        #  "orchid", "orchid1", "orchid2",
        #  "orchid3", "orchid4", "darkorchid",
        #  "darkorchid1", "darkorchid2", "darkorchid3",
        #  "darkorchid4", "purple", "purple1",
        #  "purple2", "purple3", "purple4"
        # )) {

        distance_simulated <-
          compute_cluster_splits(
            input_simulated$sample_description, input_simulated$post_sampling,
            simulation_instance, input_simulated$n_cells, input_simulated$n_mutations,
            input_simulated$n_clusters, input_simulated$allele_count,
            input_simulated$mutated_read_counts, input_simulated$total_read_counts,
            n_mutation_sampling_events = n_sampling_events,
            n_tree_sampling_events = n_sampling_events,
            cell_pair_selection = "orchid"
          )

        plot(
          ggplot(
            data.frame(x = distance_simulated$aggregatedBranchingProbabilities), aes(x = x)
          ) +
            geom_histogram(binwidth = 0.01)
        )

        mean_splitting_scores_oligo <- c(mean_splitting_scores_oligo, mean(distance_simulated$aggregatedBranchingProbabilities))



        involved_cell_indices <- sub(paste0(".*", tree_name, "_"), "", simulation_instance) %>%
          strsplit("_") %>%
          unlist()
        involved_cell_indices <- as.numeric(involved_cell_indices)

        involved_single_tumor_cells <- description_data$sample_name[involved_cell_indices]

        pairs <- combn(involved_single_tumor_cells, 2, simplify = FALSE)

        distance_separate <-
          compute_cluster_splits(input$sample_description, input$post_sampling,
            tree_name, input$n_cells, input$n_mutations,
            input$n_clusters, input$allele_count,
            input$mutated_read_counts, input$total_read_counts,
            n_mutation_sampling_events = n_sampling_events,
            n_tree_sampling_events = n_sampling_events,
            cell_pair_selection = pairs
          )

        plot(
          ggplot(
            data.frame(x = distance_separate$aggregatedBranchingProbabilities), aes(x = x)
          ) +
            geom_histogram(binwidth = 0.01)
        )
        deviance_splitting_score <- c(deviance_splitting_score, -mean(distance_separate$aggregatedBranchingProbabilities) + mean(distance_simulated$aggregatedBranchingProbabilities))
      },
      error = function(e) {
        print(e)
      }
    )
  }
}


mean_splitting_scores_mono <- mean_splitting_scores_mono[!is.na(mean_splitting_scores_mono)]


save(mean_splitting_scores_oligo, file = "~/Documents/projects/CTC_backup/simulations/simulation3_new/mean_branching_probs_oligo2.RData")
save(deviance_splitting_score, file = "~/Documents/projects/CTC_backup/simulations/simulation3_new/deviance_splitting_score2.RData")
save(mean_splitting_scores_mono, file = "~/Documents/projects/CTC_backup/simulations/simulation3_new/mean_branching_probs_mono.RData")


mean_splitting_scores_mono2 <- mean_splitting_scores_mono
mean_splitting_scores_oligo2 <- mean_splitting_scores_oligo
deviance_splitting_score2 <- deviance_splitting_score

#########

mean_splitting_scores_mono <- mean_splitting_scores_mono[!is.na(mean_splitting_scores_mono)]

data.frame(splitting_probs = mean_splitting_scores_mono, Monoclonal = TRUE) %>%
  ggplot(aes(y = splitting_probs)) +
  geom_boxplot() +
  theme_minimal()


load("~/Documents/projects/CTC_backup/simulations/simulation3_new/mean_branching_probs_oligo2.RData")
load("~/Documents/projects/CTC_backup/simulations/simulation3_new/mean_branching_probs_mono.RData")
load("~/Documents/projects/CTC_backup/simulations/simulation3_new/deviance_splitting_score2.RData")


set.seed(555)

splitting_probs <- data.frame(splitting_probs = mean_splitting_scores_mono, Oligoclonal = "Monoclonal")
splitting_probs2 <- data.frame(splitting_probs = mean_splitting_scores_oligo, Oligoclonal = "Oligoclonal")
splitting_probs_single_cells <- data.frame(splitting_probs = -deviance_splitting_score + mean_splitting_scores_oligo, Oligoclonal = "Genetically distinct single cells")

splitting_probs <- rbind(splitting_probs, splitting_probs2, splitting_probs_single_cells)
splitting_probs <- splitting_probs %>% mutate(Oligoclonal = factor(Oligoclonal, levels = c("Monoclonal", "Oligoclonal", "Genetically distinct single cells")))
levels(splitting_probs$Oligoclonal) <- c("Simulated \n monoclonal", "Simulated \n oligoclonal", "Genetically \n distinct single cells")


dim(splitting_probs[splitting_probs$Oligoclonal == "Simulated \n monoclonal", ])
dim(splitting_probs[splitting_probs$Oligoclonal == "Simulated \n oligoclonal", ])
dim(splitting_probs[splitting_probs$Oligoclonal == "Genetically \n distinct single cells", ])

final_output <- sample(rownames(splitting_probs[splitting_probs$Oligoclonal == "Simulated \n monoclonal", ]), replace = FALSE, size = 50)
final_output <- c(final_output, sample(rownames(splitting_probs[splitting_probs$Oligoclonal == "Simulated \n oligoclonal", ]), replace = FALSE, size = 50))
final_output <- c(final_output, sample(rownames(splitting_probs[splitting_probs$Oligoclonal == "Genetically \n distinct single cells", ]), replace = FALSE, size = 50))

splitting_probs %>%
  filter(rownames(.) %in% final_output) %>%
  ggplot(aes(y = splitting_probs, x = Oligoclonal, group = Oligoclonal)) +
  geom_boxplot(fill = "#41B7C4") +
  ylab("Mean splitting score") +
  xlab("Clonality status of CTC cluster") +
  theme_classic() +
  theme(
    # axis.title = element_text(size = 0),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  "/Users/jgawron/Documents/projects/CTC_backup/simulations/supplementary_figure_simulation.pdf",
  width = 6, height = 4, units = "in"
)

ggsave(
  "/Users/jgawron/Documents/projects/CTC_backup/simulations/supplementary_figure_simulation.png",
  width = 6, height = 4, units = "in"
)



data.frame(y = deviance_splitting_score) %>%
  ggplot(aes(y = y)) +
  geom_boxplot() +
  theme_minimal()


data.frame(y = -deviance_splitting_score + mean_splitting_scores_oligo) %>%
  ggplot(aes(y = y)) +
  geom_boxplot() +
  theme_minimal()
