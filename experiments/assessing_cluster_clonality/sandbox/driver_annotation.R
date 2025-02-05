library(tidyverse)

tree_names <-
  c("Br11", "Br23", "Br38", "Br39", "Br57", "Br61", "Br7", "Brx50", "Pr6", "Pr9")

for (tree in tree_names) {
  data <-
    read_delim(
      paste0("~/Documents/projects/CTC_backup/input_folder/filtered/CGI/", tree, "_cgi/alterations.tsv"),
      delim = "\t"
    )
  gene_annotation <-
    read_delim(
      paste0("~/Documents/projects/CTC_backup/input_folder/", tree, "/", tree, "_variants_annotations.csv"),
      delim = ","
    )
  gene_annotation$oncogenic_pred <- FALSE
  gene_annotation$chromosome <- NA
  gene_annotation$position <- NA


  split_gene_name <- function(gene_name) {
    split_string <- strsplit(gene_name, split = "_")[[1]]
    chromosome <- split_string[1]
    chromosome <- sub("^chr", "", chromosome)
    position <- split_string[length(split_string)]
    return(list(chromosome = chromosome, position = position))
  }



  idx <- 0
  for (gene in gene_annotation$variantName) {
    idx <- idx + 1
    genomic_pos <- split_gene_name(gene)
    gene_annotation$chromosome[idx] <- paste0("chr", genomic_pos$chromosome)
    gene_annotation$position[idx] <- genomic_pos$position
    oncogenic_pred <- data %>%
      dplyr::filter(
        CHROMOSOME == genomic_pos$chromosome, POSITION == genomic_pos$position
      ) %>%
      dplyr::select(`CGI-Oncogenic Prediction`)
    if (nrow(oncogenic_pred) == 0) {
      next
    }
    if (!(is.na(oncogenic_pred))) {
      if (grepl("driver", oncogenic_pred)) {
        gene_annotation$oncogenic_pred[idx] <- TRUE
      }
    }
  }

  gene_annotation2 <- gene_annotation[, c("chromosome", "position", "REF", "ALT")]

  write_delim(
    gene_annotation2,
    file = paste0("~/Documents/projects/CTC_backup/input_folder/", tree, "/", tree, "mutations_file.tsv"),
    delim = " "
  )


  write_csv(
    gene_annotation,
    paste0("~/Documents/projects/CTC_backup/input_folder/", tree, "/", tree, "_variants_annotations_with_driver.csv")
  )
}
