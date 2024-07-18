source("../resources/functions.R")
library("optparse")

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input-file"),
  type = "character",
  default = "~/Documents/projects/CTC_backup/input_folder", help = "Path to the folder containing all input files"
)
parser <- add_option(parser, c("-n", "--name-of-tree"),
  type = "character",
  default = "Br23", help = "Name of the tree for which to simulate CTC-clusters"
)
args <- parse_args(parser, args = c("--input-file", "--name-of-tree"))




inputFolder <- dirname(args$"input-file")
treeName <- args$name_of_tree


# inputFolder <- "~/Documents/projects/CTC_backup/input_folder"
# treeName <- "Br23"


input <- load_data(inputFolder, treeName)

allClusterSizes <- input$sample_description %>%
  filter(WBC == 0 & color != "gray93") %>%
  group_by(color) %>%
  filter(n() > 1) %>%
  summarize(cluster_size = n()) %>%
  dplyr::select("cluster_size") %>%
  unique()

write_csv(allClusterSizes, file.path(inputFolder, treeName, paste(treeName, "clusterSizes.csv", sep = "_")))
