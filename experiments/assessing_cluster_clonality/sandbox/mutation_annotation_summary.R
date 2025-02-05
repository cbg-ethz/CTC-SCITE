library(tidyverse)
source("~/Documents/projects/CTC-SCITE/experiments/assessing_cluster_clonality/workflow/scripts/simulateCTCclusters.R")

input_folder <- "~/work/ctc-data/WES_experiment/tree_sampling/"
splitting_summary <- read_delim("~/work/ctc-data/WES_experiment/splitting_summaries/splittingSummary_full_with_sample_names.tsv", delim = "\t")
splitting_summary <- splitting_summary %>% dplyr::filter(`Sample Name` != "Lu2")
tree_names <- splitting_summary %>%
  dplyr::select(`Sample Name`) %>%
  unique()
splitting_summary$mutational_burden <- NA
splitting_summary$high_impact_mutations <- NA
splitting_summary$medium_impact_mutations <- NA
splitting_summary$drivers <- NA

old_tree_name <- ""
for (idx in seq_len(nrow(splitting_summary))) {
  tree_name <- splitting_summary$`Sample Name`[idx]

  if (!tree_name == old_tree_name) {
    input <- load_data(input_folder, tree_name)


    genotypes <- call_genotypes(n_tree_sampling_events = 100, input = input)

    genotypes_wide <- genotypes %>%
      dplyr::select(Mutation, Sample, Genotype) %>%
      pivot_wider(names_from = Mutation, values_from = Genotype)

    genotypes_wide$Sample <- input$sample_description$color

    variants_annotation <- read_delim(paste0("~/work/ctc-data/WES_experiment/splitting_summaries/", tree_name, "/", tree_name, "_variants_annotations_with_driver.csv"))
  }
  old_tree_name <- tree_name

  sample_color <- splitting_summary$Color[idx]

  called_genotypes <- genotypes_wide %>%
    dplyr::filter(Sample == sample_color) %>%
    dplyr::select(-Sample) %>%
    apply(2, function(col) any(col == 1))


  called_variants_annotated <- variants_annotation %>% filter(called_genotypes)


  splitting_summary$mutational_burden <-
    called_variants_annotated %>% nrow()

  splitting_summary$medium_impact_mutations[idx] <-
    called_variants_annotated %>%
    dplyr::filter(relevant == "MODERATE") %>%
    nrow()


  splitting_summary$high_impact_mutations[idx] <-
    called_variants_annotated %>%
    dplyr::filter(relevant == "HIGH") %>%
    nrow()

  splitting_summary$drivers[idx] <- sum(called_variant_annotated$soncogenic_pred)
}

write_delim(splitting_summary, file = "~/work/ctc-data/WES_experiment/splitting_summaries/splittingSummary_full_with_sample_names_annotated.tsv", delim = "\t")
