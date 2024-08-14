## Author: Johannes Gawron
## Colors and labels in the final figure may deviate through manual manipulation in Adobe Illustrator

library(MCMCprecision)
library(ggplot2)
library(tidyverse)
library(poolr)
library(ComplexHeatmap)
library(circlize)

compute_p_value <- function(k, n, r, samplingSize, primaryTumorVector) {
  if (k > n) {
    print("Warning: Number of monoclonal clusters higher than total number of clusters!")
  }
  sampleVector <- rdirichlet(samplingSize, 1 + primaryTumorVector)
  successes <- c()
  expectedMono <- c()
  foldChange <- c()
  for (i in 1:samplingSize) {
    successRate <- sum(sampleVector[i, ]^{
      r
    })
    successesTemp <- c()
    for (j in k:n) {
      successesTemp <- c(successesTemp, dbinom(j, n, successRate))
    }
    successes <- c(successes, sum(successesTemp))
    expectedMono <- c(expectedMono, successRate * n)
    foldChange <- c(foldChange, k / (successRate * n))
  }
  return(list(p_value = sum(successes) / samplingSize, expected_mono = sum(expectedMono) / samplingSize, nFoldChange = sum(foldChange) / samplingSize))
}



summary_df_filter <- readRDS("summary_df_filter_final.rds")
summary_df_filter_combined <- summary_df_filter %>%
  mutate(group = ifelse(value %in% c(100, 1000), "100-1000",
    ifelse(value == 10000, "10000", "50000")
  ))


setwd("~/Documents/projects/CTC_backup/validation_experiment/")


p_values_table <- data.frame(matrix(ncol = 0, nrow = 24))
all_clusters <- data.frame(matrix(ncol = 0, nrow = 24))
mono_clusters <- data.frame(matrix(ncol = 0, nrow = 24))
expected_values_null_table <- data.frame(matrix(ncol = 0, nrow = 24))
fold_change_table <- data.frame(matrix(ncol = 0, nrow = 24))

mouse_models <- c("140", "141", "902", "903", "904", "905", "910")

for (mouse_model in mouse_models) {
  setwd("Cluster_csv_files")
  files <- list.files(pattern = paste0("^.*_", mouse_model, ".*\\.csv$"))
  print(files)
  summary_temp <- summary_df_filter_combined[paste0(summary_df_filter_combined$basename, ".csv") %in% files, ]


  total_number_of_clusters_by_size <- summary_temp %>%
    group_by(cluster_size) %>%
    summarize(all_clusters = n())


  number_of_monoclonal_cluster_by_size <- summary_temp %>%
    filter(clonality == "mono") %>%
    group_by(cluster_size) %>%
    summarize(monoclonal_clusters = n())

  summary_clusters_by_size <- merge(total_number_of_clusters_by_size, number_of_monoclonal_cluster_by_size, by.x = "cluster_size", all.x = TRUE)
  summary_clusters_by_size <- merge(data.frame(cluster_size = 2:25), summary_clusters_by_size, by.x = "cluster_size", all.x = TRUE)
  summary_clusters_by_size[is.na(summary_clusters_by_size)] <- 0

  setwd("../Primary_tumor_csv_files")


  primary_files <- list.files(pattern = paste0("^.*", mouse_model, ".*\\.csv$"))

  print("Loading primary tumor data:")
  primaryA <- read_delim(primary_files[1], col_names = F, delim = "\t")


  colnames(primaryA)[1] <- "barcodes"
  colnames(primaryA)[3] <- "counts"

  primaryA_counts <- primaryA %>%
    filter(counts != 0)

  if (length(primary_files) > 1) {
    primaryB <- read_delim(primary_files[2], col_names = F, delim = "\t")
    colnames(primaryB)[1] <- "barcodes"
    colnames(primaryB)[3] <- "counts"

    primaryB_counts <- primaryB %>%
      filter(counts != 0)

    length(intersect(primaryA_counts$barcodes, primaryB_counts$barcodes))

    primary_counts <- rbind(primaryA_counts, primaryB_counts)
    primary_counts <- primary_counts %>%
      filter(barcodes %in% intersect(primaryA_counts$barcodes, primaryB_counts$barcodes))
  } else {
    primary_counts <- primaryA_counts
  }

  print("Done!")
  setwd("..")

  primary_counts <- primary_counts %>%
    group_by(barcodes) %>%
    summarize(total_count = sum(counts))



  print("Compute p-values for:")
  print(summary_clusters_by_size)
  p_values_for_cluster <- rep(NA, 24)
  frequency_mono_for_cluster <- rep(NA, 24)
  nfold_change_for_cluster <- rep(NA, 24)
  expected_values_null <- rep(NA, 24)
  for (i in 1:nrow(summary_clusters_by_size)) {
    if (summary_clusters_by_size[i, 2] > 0) {
      test <- compute_p_value(summary_clusters_by_size[i, 3],
        summary_clusters_by_size[i, 2],
        summary_clusters_by_size[i, 1],
        samplingSize = 10000,
        primaryTumorVector = primary_counts$total_count
      )
      p_value <- test$p_value
      fold_change <- test$nFoldChange
      expected_value_null <- test$expected_mono
      frequency <- summary_clusters_by_size[i, 3] / summary_clusters_by_size[i, 2]
    } else {
      p_value <- NA
      frequency <- NA
      expected_value_null <- NA
      fold_change <- NA
    }
    p_values_for_cluster[i] <- p_value
    frequency_mono_for_cluster[i] <- frequency
    nfold_change_for_cluster[i] <- fold_change
    expected_values_null[i] <- expected_value_null
  }

  p_values_table <- cbind(p_values_table, data.frame(p_values_for_cluster))
  names(p_values_table)[length(names(p_values_table))] <- mouse_model

  all_clusters <- cbind(all_clusters, data.frame(summary_clusters_by_size$all_clusters))
  mono_clusters <- cbind(mono_clusters, data.frame(summary_clusters_by_size$monoclonal_clusters))
  expected_values_null_table <- cbind(expected_values_null_table, data.frame(expected_values_null))

  names(all_clusters)[length(names(all_clusters))] <- mouse_model
  names(mono_clusters)[length(names(mono_clusters))] <- mouse_model
  names(expected_values_null_table)[length(names(expected_values_null_table))] <- mouse_model

  fold_change_table <- cbind(fold_change_table, data.frame(nfold_change_for_cluster))
  names(fold_change_table)[length(names(fold_change_table))] <- mouse_model
}


frequency_table <- mono_clusters / all_clusters
rownames(frequency_table) <- 2:25
rownames(p_values_table) <- 2:25
rownames(fold_change_table) <- 2:25



frequency_table[is.na(frequency_table)] <- -1
log_fold_change_table <- log(fold_change_table + 1)


frequency_table <- t(frequency_table[1:4, ])
frequency_table <- frequency_table[c(5, 6, 2, 1, 3, 4, 7), ]

p_values_table <- t(p_values_table[1:4, ])
p_values_table <- p_values_table[c(5, 6, 2, 1, 3, 4, 7), ]

fold_change_table[is.na(fold_change_table)] <- -1
log_fold_change_table[is.na(log_fold_change_table)] <- -1


col_fun <- colorRamp2(c(-1, 0, 0.25, 0.5, 0.75, 1), c("ghostwhite", "#F7FCC9", "#EDF8BC", "#C7EBB1", "#74C9BC", "#41B7C4"))




p_values_rows <- apply(p_values_table, MARGIN = 1, FUN = function(x) {
  x <- x[!is.na(x)]
  return(fisher(x)$p)
})
p_values_columns <- apply(p_values_table, MARGIN = 2, FUN = function(x) {
  x <- x[!is.na(x)]
  return(fisher(x)$p)
})
annotation_col <- data.frame(
  "combined p-values per samples" = p_values_columns
)
annotation_row <- data.frame(
  "combined p-values per cluster size" = p_values_rows
)



columns_annotation <- c()
for (value in p_values_columns) {
  if (value < 0.001) {
    columns_annotation <- c(columns_annotation, "***")
  } else if (value < 0.01) {
    columns_annotation <- c(columns_annotation, "**")
  } else if (value < 0.05) {
    columns_annotation <- c(columns_annotation, "*")
  } else {
    columns_annotation <- c(columns_annotation, "")
  }
}

row_annotation <- c()
for (value in p_values_rows) {
  if (value < 0.001) {
    row_annotation <- c(row_annotation, "***")
  } else if (value < 0.01) {
    row_annotation <- c(row_annotation, "**")
  } else if (value < 0.05) {
    row_annotation <- c(row_annotation, "*")
  } else {
    row_annotation <- c(row_annotation, "")
  }
}

ha_col <- HeatmapAnnotation(
  foo = anno_text(x = columns_annotation, location = 1, rot = 45, show_name = TRUE, gp = gpar(fontsize = 8)),
  annotation_label = c("Sample"), annotation_name_gp = gpar(fontsize = 10)
)

ha_row <- rowAnnotation(
  foo = anno_text(x = row_annotation, which = "row", show_name = TRUE, gp = gpar(fontsize = 8)),
  annotation_label = c("Cluster size"), annotation_name_gp = gpar(fontsize = 10)
)



ht <- Heatmap(frequency_table,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(p_values_table[i, j])) {
      if (p_values_table[i, j] < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 8))
      } else if (p_values_table[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 8))
      } else if (p_values_table[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 8))
      }
    } else if (is.na(p_values_table[i, j])) {
      grid.text("-", x, y, gp = gpar(fontsize = 8))
    }
  },
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  heatmap_legend_param = list(
    title = "% mono",
    at = c(0, 0.5, 1), labels = c("0", "50", "100")
  ),
  top_annotation = ha_col,
  left_annotation = ha_row
)

lgd_sig <- Legend(pch = c("*", "**", "***", "-"), type = "points", labels = c("<0.05", "<0.01", "<0.001", "NA"), legend_gp = gpar(fontsize = 8))
draw(ht, annotation_legend_list = list(lgd_sig))

##### Same only for fold-change

col_fun <- colorRamp2(c(-1, 0, 1, 3, 5, 10), c("ghostwhite", "#F7FCC9", "#EDF8BC", "#C7EBB1", "#74C9BC", "#41B7C4"))
ht <- Heatmap(log_fold_change_table,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if (!is.na(p_values_table[i, j])) {
      if (p_values_table[i, j] < 0.001) {
        grid.text("***", x, y, gp = gpar(fontsize = 8))
      } else if (p_values_table[i, j] < 0.01) {
        grid.text("**", x, y, gp = gpar(fontsize = 8))
      } else if (p_values_table[i, j] < 0.05) {
        grid.text("*", x, y, gp = gpar(fontsize = 8))
      }
    } else if (is.na(p_values_table[i, j])) {
      grid.text("-", x, y, gp = gpar(fontsize = 8))
    }
  },
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  heatmap_legend_param = list(
    title = "log fold-change",
    at = c(0, 5, 10), labels = c("0", "5", "10")
  ),
  top_annotation = ha_col,
  left_annotation = ha_row
)

lgd_sig <- Legend(pch = c("*", "**", "***", "-"), type = "points", labels = c("<0.05", "<0.01", "<0.001", "NA"), legend_gp = gpar(fontsize = 8))
draw(ht, annotation_legend_list = list(lgd_sig))
