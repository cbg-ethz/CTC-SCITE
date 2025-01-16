library(tidyverse)
data <- read_tsv("~/Downloads/splittingSummary_full_final.tsv")
data <- data %>% mutate(WBC = n_wbcs > 0)

data <- data %>% mutate(impact_mutations = high_impact_mutations + medium_impact_mutations)


View(data)
filtered_data <- data %>%
  filter(str_detect(`Sample Name`, "Br|Pr|LM2"))


fit <- glm(as.factor(Oligoclonal) ~ n_cells + `Sample Name` + WBC, data = data, family = binomial(link = "logit"))
summary(fit)
#### Not signficant


fit2 <- glm(as.factor(Oligoclonal) ~ n_cells + `Sample Name` + WBC, data = filtered_data, family = binomial(link = "logit"))
summary(fit2)
#### Not signficant


fit3 <- glm(n_wbcs ~ n_cells + `Sample Name` + Oligoclonal, data = filtered_data, family = poisson(link = "log"))
summary(fit3)

fit4 <- glm(n_wbcs ~ n_cells + `Sample Name` + Oligoclonal, data = filtered_data, family = poisson(link = "log"))
summary(fit4)
### Not significant


fit6 <- glm(high_impact_mutations ~ n_cells + `Sample Name` + Oligoclonal + WBC, data = filtered_data, family = poisson(link = "log"))
summary(fit6)
### Not significant


fit7 <- glm(as.factor(Oligoclonal) ~ n_cells + `Sample Name` + impact_mutations + Oligoclonal + WBC, data = filtered_data, family = binomial(link = "logit"))
summary(fit7)


fit5 <- glm(impact_mutations ~ n_cells + `Sample Name` + Oligoclonal + WBC, data = filtered_data, family = poisson(link = "log"))
summary(fit5)
### significant effect of WBC presence on impact mutations

filtered_data %>%
  ggplot(aes(x = WBC, y = impact_mutations, group = WBC)) +
  geom_boxplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 18)
  ) +
  ylab("# Mutations with functional impact")



fit8 <- glm(impact_mutations ~ n_cells + `Sample Name` + Oligoclonal + n_wbcs, data = filtered_data, family = poisson(link = "log"))
summary(fit8)


filtered_data %>%
  ggplot(aes(x = n_wbcs, y = impact_mutations, group = n_wbcs)) +
  geom_boxplot() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 18)
  ) +
  ylab("# Mutations with functional impact")



fit9 <- glm(n_cells ~ `Sample Name` + Oligoclonal, data = filtered_data, family = poisson(link = "log"))
summary(fit9)
### Not significant
