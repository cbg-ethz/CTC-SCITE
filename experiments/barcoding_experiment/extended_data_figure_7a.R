library(tidyverse)
library(boot)
# library(smoothr)

setwd("~/work/ctc-data/barcoding_experiment/combined_primary_cluster/")

tables <- list()

filelist <- list.files(path = "~/work/ctc-data/barcoding_experiment/combined_primary_cluster/", pattern = "^combined.*_filter_merge\\.rds$")

for (file in filelist) {
  tables[[file]] <- readRDS(file)
  tables[[file]]$observations <- paste0(file, tables[[file]]$observations)
  first_row_with_nonzero_counts <- tables[[file]][tables[[file]]$freq != 0, ][1, ]
  total_counts <- first_row_with_nonzero_counts$counts / first_row_with_nonzero_counts$freq
  tables[[file]]$total_counts <- total_counts
  cutoff <- quantile(tables[[file]]$prop_av, prob = 0.01)
  tables[[file]] <- tables[[file]] %>%
    filter(prop_av > cutoff) %>%
    filter(prop_av != 0)
}


# Merge all tables
merged_table <- do.call(rbind, tables)
unique(merged_table$total_counts)

merged_table <- merged_table %>% mutate(complexity = as.numeric(complexity))

summary(merged_table)



merged_table2 <- merged_table %>%
  group_by(observations) %>%
  slice(rep(1:n(), total_counts)) %>%
  mutate(present = as.integer(row_number() <= first(counts))) %>%
  ungroup()


color <- "#3C8181"

ggplot(merged_table2, aes(x = prop_av, y = present)) +
  geom_jitter(width = 0, height = 0.1, color = color) +
  geom_smooth(
    method = "glm",
    method.args = list(family = "binomial"), formula = y ~ logit(x), color = color
  ) +
  labs(
    x = "Clonal frequency in primary tumor",
    y = "P(barcode is present in CTC cluster)"
  ) +
  geom_abline(linetype = "dotted", slope = 1, intercept = 0)


get_mean_quantile <- function(cutoff, merged_table, upper = TRUE) {
  if (upper) {
    return(mean(merged_table$prop_av[merged_table$prop_av > cutoff]))
  } else {
    return(mean(merged_table$prop_av[merged_table$prop_av < cutoff]))
  }
}

get_mean_quantile2 <- function(cutoff, merged_table, upper = TRUE) {
  if (upper) {
    return(mean(merged_table$freq[merged_table$prop_av > cutoff]))
  } else {
    return(mean(merged_table$freq[merged_table$prop_av < cutoff]))
  }
}


boot_wrapper <- function(data, indices, cutoff, upper = TRUE) {
  # Subset the data based on bootstrap sample indices
  sampled_data <- data[indices, ]
  # Call the target function
  return(get_mean_quantile2(cutoff, sampled_data, upper))
}

# Perform the bootstrap
set.seed(123) # For reproducibility
cutoff <- merged_table$prop_av %>% quantile(0.999) # 0.1
bootstrap_replicates <- 1000
bootstrap_results <- boot(
  data = merged_table,
  statistic = function(data, indices) boot_wrapper(data, indices, cutoff, upper = TRUE),
  R = bootstrap_replicates # Number of bootstrap replicates
)
# Print bootstrap results
print(bootstrap_results)

bootstrap_statistics <- bootstrap_results$t

# Convert to a data frame for ggplot2
bootstrap_df <- data.frame(stats = bootstrap_statistics, bias = bootstrap_statistics - get_mean_quantile(cutoff, merged_table, upper = TRUE), upper = TRUE, mean_fraction_primary = get_mean_quantile(cutoff, merged_table, upper = TRUE))

bootstrap_results <- boot(
  data = merged_table,
  statistic = function(data, indices) boot_wrapper(data, indices, cutoff, upper = FALSE),
  R = bootstrap_replicates # Number of bootstrap replicates
)
# Print bootstrap results
print(bootstrap_results)

bootstrap_statistics <- bootstrap_results$t


bootstrap_df2 <- data.frame(stats = bootstrap_statistics, bias = bootstrap_statistics - get_mean_quantile(cutoff, merged_table, upper = FALSE), upper = FALSE, mean_fraction_primary = get_mean_quantile(cutoff, merged_table, upper = FALSE))
bootstrap_df <- rbind(bootstrap_df, bootstrap_df2)

bootstrap_df$mean_fraction_primary <- as.factor(bootstrap_df$mean_fraction_primary)
mean_frac_primary_upper <- levels(bootstrap_df$mean_fraction_primary)[2]

ggplot(bootstrap_df, aes(x = upper, y = stats)) +
  geom_boxplot(fill = "#41B7C4", color = "#3C8181") +
  labs(
    # title = "Shift in abundancy of highly against lowly represented clones",
    y = "Mean frequency among CTC clusters",
    x = ""
  ) +
  theme_classic() +
  theme(
    # axis.title = element_text(size = 0),
    legend.title = element_text(size = 0),
    legend.text = element_text(size = 0),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  annotate(
    "segment",
    x = 1.5,
    xend = 2.5,
    y = as.numeric(mean_frac_primary_upper),
    yend = as.numeric(mean_frac_primary_upper),
    linetype = "dashed",
    color = "red"
  ) +
  annotate(
    "text",
    x = 2,
    y = as.numeric(mean_frac_primary_upper) - 0.01,
    label = "Mean frequency in primary tumor",
    color = "red",
    size = 4,
    fontface = "italic"
  ) +
  scale_x_discrete(
    labels = c("Lowly abundant clones", "Highly abundant clones")
  )
# fill = "#41B7C4", color ="#3C8181"

ggsave(
  "~/work/ctc-data/barcoding_experiment/extended_data_figure_7a.pdf",
  width = 6, height = 4, units = "in"
)

merged_table %>%
  ggplot(aes(y = freq, x = prop_av)) +
  geom_point() +
  labs(x = "tumor VAF", y = "fraction among CTC clusters")



cor(x = merged_table2$prop_av, y = merged_table2$present, method = "spearman")
cor.test(
  x = merged_table2$prop_av, y = merged_table2$present, method = "spearman"
)
