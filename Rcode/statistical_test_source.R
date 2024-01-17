library(readr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(DescTools, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(purrr, quietly = TRUE)




reverse_paste <- function(filename, string){
  paste0(string,filename)
}


load_cluster_df <- function(filename){
  df <- read_delim(filename, delim = "\t",col_names = FALSE)
  colnames(df)[1] <- "barcodes"
  colnames(df)[3] <- "counts"
  df <- df %>% arrange(desc(counts))
  df <- df %>%
    mutate(prop_col = counts/sum(counts), cumprop_col = cumsum(prop_col))
  return(df)
}


#' Load cluster files.
#' 
#' load_cluster_files()  laods all the CTC-cluster csv-files from the validation
#' experiment and respecitve summaries intointo memory.
#' 
#'
#' @return A named list.
#' summary: A data frame which contains sumamry statistics
#' cluster_data: a list. Each entry is a data frame with 4 colums: The first 
#'  is a unique barcode identifier, and the third is the total read count.
#' 
#' @export
#'
#' @examples load_cluster_files()
#' @noRD
load_cluster_files <- function(){
  
  
  files <- list.files(path = "../../validation_experiment/Cluster_csv_files/", pattern = "\\.csv$")
  

  
  files <- map_chr(files, reverse_paste, "../../validation_experiment/Cluster_csv_files/")
  
  myfiles <- map(files, load_cluster_df)
  
  
  
  for (file in files) {
    basename <- gsub("\\.csv", "", file)
    tryCatch({
      df <- read_delim(paste0("../../validation_experiment/Cluster_csv_files/",file), delim = "\t",col_names = F)
      colnames(df)[1] <- "barcodes"
      colnames(df)[3] <- "counts"
      df_sorted <- df %>% arrange(desc(counts))
    
      #df_sorted <- df[order(-df[, 3]), ] #Sort data frame for decreasing read counts
      df_sorted <- df_sorted %>%
        mutate(prop_col = counts/sum(counts), cumprop_col = cumsum(prop_col))
      # Calculate the read proportion for each barcode and
      # the cumulative proportion
    
      df_list[[basename]] <- df_sorted
    }, error = function(e) {
      cat("Error reading ", file, "- skipping\n")
    })
  }


# create empty data frame to store summary information
summary_df <- data.frame(basename = character(),
                         num_rows_accumulating_nine = integer(),
                         num_rows_accumulating_ninefive = integer(),
                         stringsAsFactors = FALSE)

for (i in seq_along(df_list)) {
  df <- df_list[[i]]
  cumprop_col <- df[, "cumprop_col"]
  num_rows_accumulating_nine <- sum(cumprop_col <= 0.90)
  num_rows_accumulating_ninefive <- sum(cumprop_col <= 0.95)
  basename <- names(df_list)[i]
  # determine the value based on the basename
  value <- switch(substr(basename, 1, 4),
                  "10k_" = 10000,
                  "50k_" = 50000,
                  "100_" = 100,
                  "1000" = 1000,
                  NA) # if no match, assign NA
  # determine if num_rows_accumulating > 1
  more_than_one <- ifelse(num_rows_accumulating_ninefive > 0, "Yes", "No")
  prop_col_1 <- df$prop_col[1]
  prop_col_2 <- df$prop_col[2]
  summary_row <- data.frame(basename = basename,
                            num_rows_accumulating_nine = num_rows_accumulating_nine,
                            num_rows_accumulating_ninefive = num_rows_accumulating_ninefive,
                            prop_col_1 = prop_col_1,
                            prop_col_2 = prop_col_2,
                            value = value,
                            more_than_one = more_than_one,
                            stringsAsFactors = FALSE)
  summary_df <- rbind(summary_df, summary_row)
  }
  save(summary_df, file="../../validation_experiment/output/summary_df_nine_ninefive.rds")
  return(list("summary" = summary_df, "cluster_data" = df_list))
}

input_data <- load_cluster_files()


#################
####Debugging####
#################


# create empty list to store data frames
df_list <- list()
# loop through each file, load into data frame, and add to list



basename <- gsub("\\.csv", "", "50k_7_910_b_S276_R2_001.fastq.gz_stats.csv")
df <- read_delim(paste0("../../validation_experiment/Cluster_csv_files/",file), delim = "\t",col_names = F)
colnames(df)[1] <- "barcodes"
colnames(df)[3] <- "counts"
df_sorted <- df %>% arrange(desc(counts))
df_sorted <- df_sorted %>%
  mutate(prop_col = counts/sum(counts), cumprop_col = cumsum(prop_col))
df_list[[basename]] <- df_sorted

for (file in files) {
  basename <- gsub("\\.csv", "", file)
  tryCatch({
    df <- read_delim(paste0("../../validation_experiment/Cluster_csv_files/",file), delim = "\t",col_names = F)
    colnames(df)[1] <- "barcodes"
    colnames(df)[3] <- "counts"
    df_sorted <- df %>% arrange(desc(counts))
    
    #df_sorted <- df[order(-df[, 3]), ] #Sort data frame for decreasing read counts
    df_sorted <- df_sorted %>%
      mutate(prop_col = counts/sum(counts), cumprop_col = cumsum(prop_col))
    # Calculate the read proportion for each barcode and
    # the cumulative proportion
    
    df_list[[basename]] <- df_sorted
  }, error = function(e) {
    cat("Error reading ", file, "- skipping\n")
  })
}

##############






summary_df$value[summary_df$basename %like% "^1000"] <- 1000
summary_df$value[is.na(summary_df$value)] <- 10000

summary_df$cluster_size <-            ifelse(grepl("_0_", summary_df$basename), "0",
                                             ifelse(grepl("_1_", summary_df$basename), "1",
                                                    ifelse(grepl("_2_", summary_df$basename), "2",
                                                           ifelse(grepl("_3_", summary_df$basename), "3",
                                                                  ifelse(grepl("_4_", summary_df$basename), "4",
                                                                         ifelse(grepl("_5_", summary_df$basename), "5",
                                                                                ifelse(grepl("_6_", summary_df$basename), "6",
                                                                                       ifelse(grepl("_7_", summary_df$basename), "7",
                                                                                              ifelse(grepl("_8_", summary_df$basename), "8",
                                                                                                     ifelse(grepl("_9_", summary_df$basename), "9",
                                                                                                            ifelse(grepl("_10_", summary_df$basename), "10",
                                                                                                                   ifelse(grepl("_10plus_", summary_df$basename), "11",
                                                                                                                          ifelse(grepl("_11_", summary_df$basename), "11",
                                                                                                                                 ifelse(grepl("_12_", summary_df$basename), "12",
                                                                                                                                        ifelse(grepl("_13_", summary_df$basename), "13",
                                                                                                                                               ifelse(grepl("_14_", summary_df$basename), "14",
                                                                                                                                                      ifelse(grepl("_20_", summary_df$basename), "20",
                                                                                                                                                             ifelse(grepl("_25_", summary_df$basename), "25",
                                                                                                                                                                    NA))))))))))))))))))

summary_df$cluster_category_no_WBCs <- "None"

for (i in c(0:15,25)){
  pattern <- paste0("_",paste0(as.character(i),"_"))
  summary_df$cluster_category_no_WBCs[grepl(pattern, summary_df$basename)] <- as.character(i)
}


#summary_df$cluster_category_no_WBCs <-  ifelse(grepl("_0_", summary_df$basename), "0",
#                                              ifelse(grepl("_1_", summary_df$basename), "1",
#                                                    ifelse(grepl("_2_", summary_df$basename), "2",
#                                                          ifelse(grepl("_3_", summary_df$basename), "3",
#                                                                ifelse(grepl("_4_", summary_df$basename), "4",
#                                                                      ifelse(grepl("_5_|_6_|_7_|_8_|_9_|_10_|_11_|_12_|_13_|_14_|_20_|_25_|10plus", summary_df_filter$basename), "5+", NA))))))




duplicated_rows <- duplicated(substr(summary_df$basename, 1, 15))
summary_df_filter <- summary_df[!duplicated_rows, ]


summary_df_filter$cluster_size <- as.numeric(summary_df_filter$cluster_size)
summary_df_filter <- summary_df_filter[summary_df_filter$num_rows_accumulating_nine <= summary_df_filter$cluster_size, ]

colors <- c("#8491B499", "#3C548899", "#8491B4FF", "#3C5488FF", "#B09C8599")

prop_comp <- summary_df_filter %>%
  group_by(value) %>%
  summarize(num_more_than_one = sum(more_than_one == "Yes"), num_total = sum(more_than_one == "No"))

y <- prop_comp$num_more_than_one
z <- prop_comp$num_more_than_one + prop_comp$num_total
