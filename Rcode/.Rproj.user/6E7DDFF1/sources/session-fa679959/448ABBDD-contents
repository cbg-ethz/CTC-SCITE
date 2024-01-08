files <- list.files(path = "../../validation_experiment/Cluster_csv_files/", pattern = "\\.csv$")
# create empty list to store data frames
df_list <- list()
# loop through each file, load into data frame, and add to list

#renv::install("readr")
#renv::install("tidyr")
#renv::install("dplyr")
#renv::install("DescTools")
#renv::install("ggplot2")
library(readr)
library(tidyr)
library(dplyr)
library(DescTools)
library(ggplot2)

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
    cat("Error reading ", e, "- skipping\n")
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

#comp_1 <- prop.test(x = c(y[1], y[2]), n = c(z[1], z[2]))
#comp_2 <- prop.test(x = c(y[2], y[3]), n = c(z[2], z[3]))
#comp_3 <- prop.test(x = c(y[3], y[4]), n = c(z[3], z[4]))

#p_val <- c(comp_1$p.value, comp_2$p.value, comp_3$p.value)

# Group the rows by the "value" column only
summary_df_filter %>%
  group_by(value) %>%
  summarize(prop = mean(more_than_one == "Yes" & !is.na(more_than_one))) %>%
  ungroup() %>% 
  
  ggplot(aes(x = as.factor(value), y = prop)) + 
  geom_bar(stat = "identity", fill = colors[1:4], width = 0.8, position = position_dodge(width=0.6)) +
  labs(x = "Primary Tumor complexity", y = "Proportion of oligoclonal Clusters")+
  theme_classic()#+
  #geom_signif(comparisons = list(c("100", "1000"), c("1000", "10000"), c("10000", "50000")), y_position = c(1, 1, 1), tip_length = 0.01, textsize = 3, 
              annotations = ifelse(p_val < 0.001, "***", 
                                   ifelse(p_val < 0.01, "**", 
                                          ifelse(p_val < 0.05, "*", "ns"))))


#Check oligoclonality of Clusters for the lowest complexity by number of cells in cluster

centi_df <- summary_df_filter[summary_df_filter$value==100, ]

#centi_df %>%
#  group_by(cluster_category_no_WBCs) %>%
#  summarize(prop = mean(more_than_one == "Yes" & !is.na(more_than_one))) %>%
#  ungroup() %>%
  
#  ggplot(aes(x = reorder(cluster_category_no_WBCs, +prop), y = prop)) +
#  geom_bar(stat = "identity", position = "dodge", fill = colors[1:4]) +
#  labs(x = "Cluster size", y = "Proportion oligoclonal Clusters")+
#  theme_classic()


centi_df %>%
  group_by(cluster_category_no_WBCs) %>%
  summarize(prop = mean(more_than_one == "Yes" & !is.na(more_than_one))) %>%
  ungroup() %>%
  mutate(cluster_category_no_WBCs = as.numeric(cluster_category_no_WBCs)) %>%
  ggplot(aes(x = cluster_category_no_WBCs, y = prop)) +
  geom_line() +
  labs(x = "Cluster size", y = "Proportion oligoclonal Clusters")+
  theme_classic()


centi_df$cluster_category_no_WBCs <- factor(centi_df$cluster_category_no_WBCs, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"))


centi_df %>%
  group_by(cluster_category_no_WBCs) %>%
  summarize(prop = mean(more_than_one == "Yes" & !is.na(more_than_one))) %>%
  ungroup() %>%
  
  filter(cluster_category_no_WBCs %in% c("1","2","3","4","5","6","7","8")) %>%
  ggplot(aes(x = cluster_category_no_WBCs, y = prop, fill = cluster_category_no_WBCs)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster size", y = "Proportion oligoclonal Clusters")+
  theme_classic()

### Now take the fraction of each of the clones in the tumor and compute the 
# theoretical probability to obtain an oligoclonal CTC-cluster to 
# under a well-mixed model

primary_904A <- read_delim("../../validation_experiment/Primary_tumor_csv_files/20220812.B-904A_R2_stats.csv",
                           col_names = F, delim = "\t")
primary_904B <- read_delim("../../validation_experiment/Primary_tumor_csv_files/20220812.B-904B_R2_stats.csv",
                           col_names = F, delim = "\t" )

colnames(primary_904A)[1] <- "barcodes"
colnames(primary_904A)[3] <- "counts"

colnames(primary_904B)[1] <- "barcodes"
colnames(primary_904B)[3] <- "counts"


primary_904A_counts <-  primary_904A %>%
  filter(counts !=0)
primary_904B_counts <-  primary_904B %>%
  filter(counts !=0) 
 
length(intersect(primary_904A_counts$barcodes,primary_904B_counts$barcodes))

primary_904_counts <- rbind(primary_904A_counts,primary_904B_counts)
primary_904_counts <- primary_904_counts %>%
  filter(barcodes %in% intersect(primary_904A_counts$barcodes,primary_904B_counts$barcodes))

primary_904_counts <- primary_904_counts %>%
  group_by(barcodes) %>%
  summarize(total_count = sum(counts))

primary_904_counts<- primary_904_counts %>%
  mutate(relative_count = total_count/sum(total_count))


primary_904_counts$relative_count_squared <- primary_904_counts$relative_count^2

probability_of_monoclonality <- c()
for (i in 2:9){
  probability_of_monoclonality <- c(probability_of_monoclonality, sum(primary_904_counts$relative_count^i))
}

theoretical_result <- data.frame(cluster_size = 2:9,
                                 monoclonality = probability_of_monoclonality,
                                 oligoclonality = 1-probability_of_monoclonality)
  
ggplot(theoretical_result,aes(x = cluster_size, y = oligoclonality)) +
  geom_line()
 

centi_df %>%
  filter(grepl("_904_", centi_df$basename)) %>%
  group_by(cluster_category_no_WBCs) %>%
  filter(n() > 1) %>%
  summarize(prop = mean(more_than_one == "Yes" & !is.na(more_than_one))) %>%
  ungroup() %>%
  mutate(cluster_category_no_WBCs = as.numeric(cluster_category_no_WBCs)) %>%
  ggplot(aes(x = cluster_category_no_WBCs, y = prop)) +
  geom_point() +
  labs(x = "Cluster size", y = "Proportion oligoclonal Clusters")+
  theme_classic() +
  geom_line(theoretical_result[1:5,], mapping = aes(x = cluster_size, y = oligoclonality))

ggsave("../../validation_experiment/output/observed_oligoclonality_vs_predicted_904.png", plot = last_plot(), dpi = 300)




###### Do the same with the next tumor
primary_905A <- read_delim("../../validation_experiment/Primary_tumor_csv_files/20220812.B-905A_R2_stats.csv",
                           col_names = F, delim = "\t")
primary_905B <- read_delim("../../validation_experiment/Primary_tumor_csv_files/20220812.B-905B_R2_stats.csv",
                           col_names = F, delim = "\t" )

colnames(primary_905A)[1] <- "barcodes"
colnames(primary_905A)[3] <- "counts"

colnames(primary_905B)[1] <- "barcodes"
colnames(primary_905B)[3] <- "counts"


primary_905A_counts <-  primary_905A %>%
  filter(counts !=0)
primary_905B_counts <-  primary_905B %>%
  filter(counts !=0) 

length(intersect(primary_905A_counts$barcodes,primary_905B_counts$barcodes))

primary_905_counts <- rbind(primary_905A_counts,primary_905B_counts)
primary_905_counts <- primary_905_counts %>%
  filter(barcodes %in% intersect(primary_905A_counts$barcodes,primary_905B_counts$barcodes))

primary_905_counts <- primary_905_counts %>%
  group_by(barcodes) %>%
  summarize(total_count = sum(counts))

primary_905_counts<- primary_905_counts %>%
  mutate(relative_count = total_count/sum(total_count))

summary(primary_905_counts)


theoretical_result_905 <- data.frame(cluster_size = 2:9,
                                 monoclonality = probability_of_monoclonality,
                                 oligoclonality = 1-probability_of_monoclonality)


centi_df %>%
  filter(grepl("_905_", centi_df$basename)) %>%
  group_by(cluster_category_no_WBCs) %>%
#  filter(n() > 1) %>%
  summarize(prop = mean(more_than_one == "Yes" & !is.na(more_than_one))) %>%
  ungroup() %>%
  mutate(cluster_category_no_WBCs = as.numeric(cluster_category_no_WBCs)) %>%
  ggplot(aes(x = cluster_category_no_WBCs, y = prop)) +
  geom_point() +
  labs(x = "Cluster size", y = "Proportion oligoclonal Clusters")+
  theme_classic() +
  geom_line(theoretical_result_905, mapping = aes(x = cluster_size, y = oligoclonality))




ggplot(theoretical_result,aes(x = cluster_size, y = oligoclonality)) +
  geom_line()






#save(summary_df, file="summary_df_nine_ninefive.rds")



