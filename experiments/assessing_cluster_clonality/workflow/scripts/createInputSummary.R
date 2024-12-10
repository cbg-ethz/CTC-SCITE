source("../resources/functions.R")
library("optparse")

parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input-folder"),
  type = "character",
  default = "~/Documents/projects/CTC_backup/input_folder", help = "Path to the folder containing all input files"
)
parser <- add_option(parser, c("-n", "--name-of-tree"),
  type = "character",
  default = "Br23", help = "Name of the tree for which to simulate CTC-clusters"
)
args <- parse_args(parser)




input_folder <- args$"input-folder"
tree_name <- args$"name-of-tree"



# input_folder <- "~/Documents/projects/CTC_backup/input_folder"
# tree_name <- "Br23"


input <- load_data(input_folder, tree_name)

all_cluster_sizes <- input$sample_description %>%
  filter(WBC == 0 & color != "gray93") %>%
  group_by(color) %>%
  filter(n() > 1) %>%
  summarize(cluster_size = n()) %>%
  dplyr::select("cluster_size") %>%
  unique()

write_csv(all_cluster_sizes, file.path(input_folder, tree_name, paste(tree_name, "clusterSizes.csv", sep = "_")))
print("Success.")