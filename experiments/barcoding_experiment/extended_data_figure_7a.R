library(ggplot2)
library(boot)
library(smoothr)

setwd("~/Downloads")

tables <- list()
filelist <- c("combined_140_order_filter_merge.rds", "combined_141_order_filter_merge.rds", "combined_902_order_filter_merge.rds", "combined_903_order_filter_merge.rds", "combined_904_order_filter_merge.rds", "combined_905_order_filter_merge.rds", "combined_910_order_filter_merge.rds")
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


fit <-
  glm(
    present ~ log(prop_av),
    data = merged_table2, family = binomial(link = "logit")
  )
fit1 <-
  glm(
    present ~ prop_av,
    data = merged_table2, family = binomial(link = "logit")
  )
fit2 <-
  glm(
    present ~ log(prop_av),
    data = merged_table2, family = binomial(link = "log")
  )
fit3 <-
  glm(present ~ logit(prop_av),
    data = merged_table2, family = binomial(link = "logit")
  )
fit4 <-
  glm(present ~ logit(prop_av),
    data = merged_table2, family = binomial(link = "log")
  )

AIC(fit)
AIC(fit1)
AIC(fit2)
AIC(fit3)
AIC(fit4)

### fit3 and fit4 have approximately the same AIC, but fit 3 transforms the x
## coordinates to a well defined range, so I will use fit3.
fit_linear_model <- lm(freq ~ log(prop_av), data = merged_table)
AIC(fit_linear_model)

merged_table %>%
  ggplot(aes(y = freq, x = prop_av)) +
  geom_point() +
  labs(x = "tumor VAF", y = "fraction among CTC clusters") +
  geom_smooth(method = "lm", formula = y ~ x)

merged_table %>%
  ggplot(aes(y = freq, x = log(prop_av))) +
  geom_point() +
  labs(x = "tumor VAF", y = "fraction among CTC clusters") +
  geom_smooth(method = "lm")


merged_table %>%
  ggplot(aes(y = logit(freq), x = logit(prop_av))) +
  geom_point() +
  labs(x = "logit(tumor VAF)", y = "logit(fraction among CTC clusters)") +
  geom_smooth(method = "lm", formula = y ~ x)



fit_logit_logit_linear <-
  lm(asin(sqrt(freq)) ~ asin(sqrt(prop_av)), data = merged_table)
plot(fit_logit_logit_linear)
primary_VAF <- seq(0.01, 0.55, 0.01)

predicted <- predict(fit4,
  newdata = data.frame(prop_av = primary_VAF),
  type = "response"
)

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


merged_table %>%
  ggplot(aes(y = freq, x = prop_av)) +
  geom_point() +
  labs(x = "tumor VAF", y = "fraction among CTC clusters") +
  geom_line(data = data.frame(prop_av = primary_VAF, freq = predicted))

means <- rep(NA, 100)
for (i in 1:25) {
  window_start <- (i - 1) / 50
  window_end <- i / 50
  means[i] <- merged_table2 %>%
    dplyr::filter(prop_av > window_start) %>%
    dplyr::filter(prop_av <= window_end) %>%
    pull(present) %>%
    mean()
}

ggplot(data.frame(y = means, x = 0:99), aes(x = x, y = y)) +
  geom_point()

ggsave(
  "~/Desktop/clonal_prevalence2.pdf",
  width = 3.6, height = 3.6, units = "in"
)

cor(x = merged_table2$prop_av, y = merged_table2$present, method = "spearman")
cor.test(
  x = merged_table2$prop_av, y = merged_table2$present, method = "spearman"
)
