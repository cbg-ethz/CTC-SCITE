# Inferring the clonality of barcoded CTC clusters
## Author:David Gremmelspacher
## Colors and labels in the final figure may deviate through manual manipulation in Adobe Illustrator



# Load libraries
library(data.table)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(webr)
library(DescTools)

# Generate a list of barcode counts files
files <- list.files(pattern = "\\.csv$")

# Create empty list to store data frames
df_list <- list()

# Loop through all sample barcode counts files
for (file in files) {
  sample_alias <- gsub("\\.csv", "", file)
  tryCatch(
    {
      # Load counts into data frame
      df <- read.delim(file)
      # Sort by decreasing barcode counts
      df_sorted <- df[order(-df[, 3]), ]
      # Calculate fraction of total counts for each Barcode
      prop_col <- df_sorted[, 3] / sum(df_sorted[, 3])
      # Calculate cumulative barcode proportions
      cumprop_col <- cumsum(prop_col)
      # Bind columns
      df_sorted <- cbind(df_sorted, prop_col, cumprop_col)
      # Add to list
      df_list[[sample_alias]] <- df_sorted
    },
    error = function(e) {
      cat("Error reading", file, "- skipping\n")
    }
  )
}

# Create empty data frame to store summary information for all samples
summary_df <- data.frame( # Sample alias
  sample_alias = character(),
  # Number of barcodes accumulating 90% of total aligned counts
  num_rows_accumulating_nine = integer(),
  # Proportion of total counts for most abundant barcode
  prop_col_1 = numeric(),
  # Proportion of total counts for second most abundant barcode
  prop_col_2 = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each sample data frame in df_list and calculate summary information to populate summmary_df
for (i in seq_along(df_list)) {
  df <- df_list[[i]]
  cumprop_col <- df[, "cumprop_col"]
  num_rows_accumulating_nine <- sum(cumprop_col < 0.9) + 1
  sample_alias <- names(df_list)[i]
  prop_col_1 <- df$prop_col[1]
  prop_col_2 <- df$prop_col[2]
  summary_row <- data.frame(
    sample_alias = sample_alias,
    num_rows_accumulating_nine = num_rows_accumulating_nine,
    prop_col_1 = prop_col_1,
    prop_col_2 = prop_col_2,
    stringsAsFactors = FALSE
  )
  summary_df <- rbind(summary_df, summary_row)
}

# Read-in metadata containing information on CTC cell count per cluster, primary tumor barcode diversity and mouse/tumor IDs
metadata <- read_excel("metadata.xlsx")

# Merge summary dataframe and metadata by CTC cluster sample_alias
merged_df <- merge(summary_df, metadata, by = "sample_alias", all = TRUE)

# Assign CTC clusters into categories "0" (negative controls), "2" or "3+" based on the CTC counts
merged_df <- merged_df %>%
  mutate(cluster_category = case_when(
    n_ctc == 0 ~ "0",
    n_ctc == 2 ~ "2",
    TRUE ~ "3+"
  ))

# Quality filtering of CTC cluster samples
merged_df$n_ctc <- as.numeric(merged_df$n_ctc)
merged_df_filter <- merged_df[merged_df$num_rows_accumulating_nine <= merged_df$n_ctc, ]

# Assign CTC cluster samples mono- or oligoclonal based on the dominance of the most abundant barcode, taking into account the cell number
merged_df_filter$clonality <- ifelse(merged_df_filter$prop_col_2 / merged_df_filter$prop_col_1 < 1 / merged_df_filter$n_ctc,
  "mono",
  "oligo"
)

# Return fraction of oligoclonal CTC clusters across all samples
nrow(merged_df_filter[merged_df_filter$clonality == "oligo", ]) / nrow(merged_df_filter)

# Assign CTC cluster samples a complexity value of "Low", "Medium" or "High", based on corresponding primary tumor barcode complexity
merged_df_filter <- merged_df_filter %>%
  mutate(complexity = ifelse(barcode_initial_complexity %in% c(100, 1000), "Low",
    ifelse(barcode_initial_complexity == 10000, "Medium", "High")
  ))


# Perform Cochran-Armitage test for trend
# Define matrix with counts for oligo- and monoclonal CTC clusters across complexities (as inferred from merged_df_filter)
x <- matrix(c(16, 133, 57, 88, 90, 42), byrow = TRUE, ncol = 2)
CochranArmitageTest(x, alternative = "one.sided")


# Fig. 2b
# Generate a data.frame with counts for mono- and oligoclonal CTC clusters within each complexity level using PieDonut
low <- data.frame(Clonality = c("Monoclonal", "Oligoclonal"), n = c(133, 16))
medium <- data.frame(Clonality = c("Monoclonal", "Oligoclonal"), n = c(88, 57))
high <- data.frame(Clonality = c("Monoclonal", "Oligoclonal"), n = c(42, 90))


# Generate a donut plot illustrating mono- and oligoclonal CTC cluster counts for each complexity level
PieDonut(low, aes(Clonality, count = n), showPieName = FALSE, donutLabelSize = 0, labelpositionThreshold = 0.01, explode = 1, r0 = 0.5, r1 = 0.95, pieLabelSize = 0, showRatioDonut = FALSE)
PieDonut(medium, aes(Clonality, count = n), showPieName = FALSE, donutLabelSize = 0, labelpositionThreshold = 0.01, explode = 1, r0 = 0.5, r1 = 0.95, pieLabelSize = 0, showRatioDonut = FALSE)
PieDonut(high, aes(Clonality, count = n), showPieName = FALSE, donutLabelSize = 0, labelpositionThreshold = 0.01, explode = 1, r0 = 0.5, r1 = 0.95, pieLabelSize = 0, showRatioDonut = FALSE)


# Create a data frame specifying for each tumor sample the proportion of oligoclonal CTC clusters in categories "2" and "3+"
combined_summary <- merged_df_filter %>%
  group_by(tumor_id, cluster_category, complexity) %>%
  summarise(
    n = n(),
    prop_oligo = mean(clonality == "oligo" & !is.na(clonality)),
    .groups = "drop"
  ) %>%
  ungroup()

# Only keep samples with counts in both categories ("2" and "3+")
combined_summary <- combined_summary[-c(21), ]

# Calculate the absolute number of mono- and oligoconal CTC clusters per category
combined_summary$oligo <- combined_summary$n * combined_summary$prop_oligo
combined_summary$mono <- combined_summary$n - combined_summary$oligo

# Sum oligoclonal CTC clusters and monoclonal CTC clusters for each category
category_2 <- subset(combined_summary, cluster_category == "2", select = c("oligo", "mono"))
category_3plus <- subset(combined_summary, cluster_category == "3+", select = c("oligo", "mono"))
total_2 <- colSums(category_2)
total_3plus <- colSums(category_3plus)

# Create the contingency table
contingency_table <- matrix(
  c(total_2["oligo"], total_2["mono"], total_3plus["oligo"], total_3plus["mono"]),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(c("Category 2", "Category 3+"), c("Oligoclonal", "Monoclonal"))
)

# Perform Fisher's Exact Test for "2" vs. "3+"
fisher.test(contingency_table)

# Generate plot for Fig. 2c
complexity_colors <- c("Low" = "#9fc8c8", "Medium" = "#54a1a1", "High" = "#1f6f6f")
ggplot(combined_summary, aes(x = as.character(cluster_category), y = prop_oligo, color = complexity)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", aes(group = tumor_id), size = 1, se = FALSE, alpha = 0.6) +
  theme_classic() +
  theme(axis.title = element_text(size = 0), legend.title = element_text(size = 0), legend.text = element_text(size = 0), axis.text = element_text(size = 0)) + # Populate labels in Adobe Illustrator
  scale_color_manual(values = complexity_colors) +
  ylim(0, 1)