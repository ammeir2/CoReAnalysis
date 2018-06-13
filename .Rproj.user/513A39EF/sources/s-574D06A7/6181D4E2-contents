library(ggplot2)
library(reshape2)
library(CoReAnalysis)

set.seed(1)
# Generating data ------
n <- 100
ncp <- rnorm(n)
df <- runif(n, 1, 120)
samp <- rt(n, df, ncp)

# Computing naive MLE ----
naive <- mapply(t_ncp_mle, samp, df)
plotdat <- data.frame(samp = samp, df = df, ncp = ncp, naive = naive)
ggplot(plotdat) + geom_point(aes(x = samp, y = naive, col = df)) +
  geom_abline(slope = 1, intercept = 1) +
  theme_bw()

# Sampling with selection ------
n <- 1000
ncp <- rnorm(n)
df <- runif(n, 10, 120)
samp <- numeric(n)
pvalThreshold <- 0.05
for(i in 1:n) {
  pval <- 1
  while(pval > pvalThreshold) {
    x <- rt(1, df = df[i], ncp = ncp[i])
    pval <- 2 * pt(-abs(x), df = df[i])
  }
  samp[i] <- x
}

# Computing conditional and naive estimates --------
naive <- mapply(t_ncp_mle, samp, df)
conditional <- mapply(t_conditional_ncp_mle, samp, df, 0.05)

plotdat <- data.frame(samp = samp, df = df, ncp = ncp,
                      naive = naive, conditional = conditional)
plotdat <- melt(plotdat, id = c("samp", "df", "ncp"))
names(plotdat)[4:5] <- c("method", "estimate")
ggplot(plotdat) +
  geom_point(aes(x = samp, y = estimate, col = df)) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  facet_wrap(~ method)

# Computing CIs ----
naive <- mapply(t_conditional_ncp_ci, x = samp, df = df,
                confidence_level = 0.95, threshold = 10^-5)
conditional <- mapply(t_conditional_ncp_ci, x = samp, df = df,
                      confidence_level = 0.95, pval_threshold = 0.05)

naive_cover <- numeric(n)
conditional_cover <- numeric(n)
for(i in 1:n) {
  naive_cover[i] <- naive[1, i] < ncp[i] & naive[2, i] > ncp[i]
  conditional_cover[i] <- conditional[1, i] < ncp[i] & conditional[2, i] > ncp[i]
}

naivedat <- data.frame(samp = samp, df = df, ncp = ncp,
                       cover = naive_cover,
                       lci = naive[1, ],
                       uci = naive[2, ],
                       method = "naive")
conddat <- data.frame(samp = samp, df = df, ncp = ncp,
                      cover = conditional_cover,
                      lci = conditional[1, ],
                      uci = conditional[2, ],
                      method = "conditional")
plotdat <- rbind(conddat, naivedat)

ggplot(plotdat) +
  geom_segment(aes(y = lci, yend = uci, x = samp, xend = samp, col = factor(cover))) +
  # geom_point(aes(x = samp, y = ncp, fill = df), size = 0.5, shape = 21) +
  theme_bw() +
  facet_wrap(~ method) +
  ylab("CI") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)

ggplot(plotdat) +
  geom_segment(aes(y = lci - ncp, yend = uci - ncp, x = samp, xend = samp, col = factor(cover))) +
  # geom_point(aes(x = samp, y = ncp, fill = df), size = 0.5, shape = 21) +
  theme_bw() +
  facet_wrap(~ method) +
  geom_hline(yintercept = 0) +
  ylab("CI - true_ncp")

# Testing wrapper and power functions ---
obs <- 3
x <- samp[obs]
obs_df <- df[obs]
pval_threshold <- 0.05
fit <- CoReAnalysis(x, test_statistic = "t", pval_threshold = 0.05,
                    df1 = obs_df, confidence_level = 1 - 0.5)
replication_power(fit, plot = TRUE)
replication_power(fit, plot = FALSE)


