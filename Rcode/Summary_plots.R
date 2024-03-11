data <- data.frame(sample = c('Br11','Br11','Br11','Br23','Br23','Br23','Br38','Br38','Br38','Br39','Br39','Br39','Br61','Br61','Br61', 'Brx50','Brx50','Brx50', 'LM2','LM2','LM2','Lu2','Lu2','Lu2',
                              'Pr6','Pr6','Pr6','Pr9','Pr9','Pr9'),
                   Category = c('Oligoclonal', 'Monoclonal', 'Likely oligoclonal','Oligoclonal', 'Monoclonal', 'Likely oligoclonal','Oligoclonal', 'Monoclonal', 'Likely oligoclonal','Oligoclonal', 'Monoclonal', 'Likely oligoclonal',
                                'Oligoclonal', 'Monoclonal', 'Likely oligoclonal','Oligoclonal', 'Monoclonal', 'Likely oligoclonal','Oligoclonal', 'Monoclonal', 'Likely oligoclonal','Oligoclonal', 'Monoclonal', 'Likely oligoclonal',
                                'Oligoclonal', 'Monoclonal', 'Likely oligoclonal','Oligoclonal', 'Monoclonal', 'Likely oligoclonal'),
                   Counts = c(0,0,1,0,4,0,0,0,1,0,0,1, 3,3,2,0,2,0,2,7,5,0,1,0,0,1,0,0,0,2))



library(ggthemes)

cbPalette <- c("#56B4E9", "#D55E00", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#999999", "#E69F00")


data %>%
  filter(Counts != 0) %>%
  ggplot(aes(x = "", y = Counts, fill = Category)) +
  geom_bar(stat = "identity") +
#  coord_polar("y", start = 0) +
  facet_wrap(~ sample, nrow = 2) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size =20)) +
  geom_text(aes(label = Counts), position = position_stack(vjust = 0.5), size = 7) +
  scale_fill_manual(values=cbPalette)


data %>% group_by(Category) %>% dplyr::summarize(TotalCounts = sum(Counts)) %>%
  ggplot(aes(x = "", y = TotalCounts, fill = Category)) +
  geom_bar(stat = "identity") +
#  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(size =20)) +
  geom_text(aes(label = TotalCounts), position = position_stack(vjust = 0.5), size = 10) +
  scale_fill_manual(values=cbPalette)


summary_data <- data.frame(Sample = c('Br7', 'Br11', 'Br23', 'Br26', 'Br38', 'Br39', 'Br61', 'LM2', 'Lu2', 'Pr6', 'Pr9', 'Br16_B', 'Br16_C', 'Br16_AC'), Oligoclonal_moderate_functional_impact = c(0,0,0,0,0,0,1,1,0,0,0,1,0,1), Oligoclonal_high_functional_impact = c(0,0,0,0,0,0,3,3,0,0,1,5,0,0), Likely_oligoclonal = c(3,1,2,1,1,0,3,7,0,0,1,0,0,1), No_oligoclonality_detected = c(0,0,0,0,0,1,1,3,1,1,0,20,8,20))

save.image(file = '~/Documents/CTC_backup/topSeparatingMutations/result_summary.RData')

library(tidyverse)

summary_data_long <- pivot_longer(summary_data, cols = c('Oligoclonal_moderate_functional_impact', 'Oligoclonal_high_functional_impact', 'Likely_oligoclonal', 'No_oligoclonality_detected'), names_to = 'Status')

summary_data_long %>% filter(value != 0) %>%
  ggplot(aes(x = '', y = value, fill = Status)) +
  geom_bar(stat = "identity") +
  #  coord_polar("y", start = 0) +
  facet_wrap(~ Sample, nrow = 2) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size =20)) +
  geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(values=cbPalette)


summary_data_long %>% group_by(Status) %>%
  summarize(totalCounts = sum(value)) %>%
  filter(totalCounts != 0) %>%
  ggplot(aes(x = '', y = totalCounts, fill = Status)) +
  geom_bar(stat = "identity") +
  #  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        strip.text.x = element_text(size =20)) +
  geom_text(aes(label = totalCounts), position = position_stack(vjust = 0.5), size = 7) +
  scale_fill_manual(values=cbPalette)


