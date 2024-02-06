library(readr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(DescTools, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(purrr, quietly = TRUE)

reverse_paste <- function(filename, string){
  paste0(string,filename)
}

reverse_paste <- function(filename, string){
  paste0(string,filename)
}


load_cluster_df <- function(filename){
  df <- read_delim(filename, delim = "\t",col_names = FALSE, col_select = c(1,3))
  colnames(df)[1] <- "barcodes"
  colnames(df)[2] <- "counts"
  df <- df %>% arrange(desc(counts)) %>%
    mutate(prop_col = counts/sum(counts), cumprop_col = cumsum(prop_col))
  return(df)
}


files <- list.files(path = "../../validation_experiment/Cluster_csv_files/", pattern = "\\.csv$")

#files <- map_chr(files, reverse_paste, "../../validation_experiment/Cluster_csv_files/")

myfiles <- map(files, load_cluster_df)


names(myfiles)<- files
save(myfiles, file = "../../validation_experiment/validation_cluster_data.Rdata")
load("../../validation_experiment/validation_cluster_data.Rdata")
names(myfiles)
