library(ggplot2)
library(CoReAnalysis)
library(dplyr)
library(magrittr)
# F and chi-square distributions
df <- 5
range <- seq(from = 0.1, to = 0.9, by = 0.1)
x <- qchisq(range, df = df, ncp = 3)
y <- qf(range, df1 = df, df2 = Inf, ncp = 3) * df
rbind(x, y)

# Generating selected data
n <- 1000
df1 <- sample(c(2, 5, 10), n, TRUE)
df2 <- sample(c(10, 50, 200), n, TRUE)
pval_threshold <- 0.05
ncp <- pmax(rf(n, df1, df2), 0)
samp <- numeric(n)
for(i in 1:n) {
  pval <- 1
  while(pval > pval_threshold) {
    x <- rf(1, df1[i], df2[i], ncp = ncp[i])
    pval <- 1 - pf(x, df1[i], df2[i])
  }
  samp[i] <- x
}

# Estimating
naive <- mapply(f_conditional_ncp_mle, x = samp, df1 = df1, df2 = df2, threshold = 0)
conditional <- mapply(f_conditional_ncp_mle, x = samp, df1 = df1, df2 = df2,
                      pval_threshold = pval_threshold)
naive <- data.frame(x = samp, df1 = df1, df2 = df2, pval_threshold = pval_threshold,
                    estimate = naive, method = "naive")
conditional <- data.frame(x = samp, df1 = df1, df2 = df2, pval_threshold = pval_threshold,
                          estimate = conditional, method = "conditional")
forplot <- rbind(naive, conditional)
ggplot(forplot, aes(x = x, y = estimate, col = method)) +
  geom_point() +
  facet_grid(df1 ~ df2, labeller = "label_both", scales = "free") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)

# Confidence intervals
naive <- mapply(f_conditional_ncp_ci, x = samp, df1 = df1, df2 = df2,
                threshold = 0, confidence_level = 0.9)
conditional <- mapply(f_conditional_ncp_ci, x = samp, df1 = df1, df2 = df2,
                      pval_threshold = 0.05, confidence_level = 0.9)
naive_cover <- numeric(n)
conditional_cover <- numeric(n)
for(i in 1:n) {
  naive_cover[i] <- naive[1, i] < ncp[i] & naive[2, i] > ncp[i]
  conditional_cover[i] <- conditional[1, i] < ncp[i] & conditional[2, i] > ncp[i]
}
mean(naive_cover)
mean(conditional_cover)

naivedat <- data.frame(samp = samp,
                       df1 = df1, df2 = df2,
                       ncp = ncp,
                       cover = naive_cover,
                       lci = naive[1, ],
                       uci = naive[2, ],
                       method = "naive")
conddat <- data.frame(samp = samp,
                      df1 = df1, df2 = df2,
                      ncp = ncp,
                      cover = conditional_cover,
                      lci = conditional[1, ],
                      uci = conditional[2, ],
                      method = "conditional")
plotdat <- rbind(conddat, naivedat)

ggplot(plotdat) +
  geom_segment(aes(y = lci, yend = uci, x = samp, xend = samp,
                   col = factor(cover))) +
  theme_bw() +
  facet_grid(method ~ df1 + df2, labeller = "label_both", scales = "free") +
  ylab("CI") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)

# Wrapper function -----
obs <- 11
x <- samp[obs]
obs_df1 <- df1[obs]
obs_df2 <- df2[obs]
pval_threshold <- 0.05
fit <- CoReAnalysis(x, test_statistic = "f", pval_threshold = 0.05,
                    df1 = obs_df1, df2 = obs_df2,
                    confidence_level = 0.95)
replication_power(fit, plot = TRUE)
replication_power(fit, plot = FALSE)
fit$threshold

