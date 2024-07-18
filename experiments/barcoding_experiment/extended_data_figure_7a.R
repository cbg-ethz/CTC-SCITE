library(tidyverse)

data <- readRDS("combi_df.rds")
View(data)


data <- data %>% mutate(counts / freq, complexity = as.numeric(complexity))

summary(data)

data2 <- data %>%
  group_by(observations) %>%
  slice(rep(1:n(), counts / freq)) %>%
  mutate(present = row_number() <= first(counts)) %>%
  ungroup()

fit <- glm(present ~ log(prop_av), data = data2, family = binomial(link = "log"))

summary(fit)

confint(fit)
coef(fit)

data %>%
  ggplot(aes(y = log(freq), x = log(prop_av))) +
  geom_point() +
  labs(x = "log(tumorVAF)", y = "log(fraction among CTC clusters)") +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2])


data %>%
  ggplot(aes(y = log(ratio), x = log(prop_av))) +
  geom_point() +
  geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2] - 1) +
  labs(x = "log(VAF in primary tumor)", y = "log(fold-change in prevalence)") +
  theme_minimal()




pdf("~/Desktop/clonal_prevalence.pdf", width = "3.6", height = "3.6")
par(bty = "l")
plot(
  x = log(data$prop_av),
  y = log(data$ratio),
  cex = 0.8,
  pch = 19,
  col = "#81A9A9",
  xlab = "log( VAF in primary tumor )",
  ylab = "log( fold-change in prevalence in CTC clusters )"
)
abline(b = coef(fit)[2] - 1, a = coef(fit)[1], col = "#3C8181", lwd = 2.2)
text(x = max(log(data$prop_av)) - 5, y = max(log(data$ratio)) - 3.5, labels = sprintf("slope = %f***", coef(fit)[2] - 1), pos = 4, col = "black")
dev.off()

summary(fit)

#### The p-value is computed manually. Compare (Regression, Modelle, Methoden
## und Anwendungen, second edition, page 204f; https://www.uni-goettingen.de/de/551357.html)

### The coefficient scaled by the standard error and then squared is
## approximately xi-square (1) distributed with one degree of freedom.

### Also note that the number of data points 1798 is large enough (>50) for a
## xi-square approximation to work, according to rule of thumb.

pseudo_t_value <- as.numeric(coef(fit)["log(prop_av)"] - 1) / as.numeric(summary(fit)$coefficients["log(prop_av)", "Std. Error"])

wald_statistics <- pseudo_t_value^2 ### xi^2_1 distributed

p.value <- 1 - pchisq(wald_statistics, 1)
print(p.value)
cor(log(data$prop_av), log(data$ratio), method = "pearson")
### The p-value is smaller that measurable by machine precision, so I suggest to report
### the smallest p-value that the generic glm summary can output:

### p-value < 2.2e-16
