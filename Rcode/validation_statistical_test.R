source("statistical_test_source.R")

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



