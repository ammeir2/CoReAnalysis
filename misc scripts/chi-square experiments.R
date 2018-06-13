library(magrittr)
library(ggplot2)

# Generating data ------
n <- 1000
df <- sample(1:100, n, replace = TRUE)
ncp <- pmax(0, rchisq(n, df = df) - df)
samp <- rchisq(n, df, ncp)

# Computing the naive MLE ------
mle <- mapply(chisq_ncp_mle, x = samp, df = df, normalize = FALSE)
plot(samp, mle)

# Sampling conditional sample -----
n <- 1000
df <- sample(c(2, 5, 20), n, replace = TRUE)
ncp <- pmax(0, rchisq(n, df = df) - df)
pval_threshold <- 0.05
samp <- numeric(n)
pvals <- numeric(n)
for(i in 1:n) {
  pvals[i] <- 1
  while(pvals[i] > pval_threshold) {
    x <- rchisq(1, df = df[i])
    pvals[i] <- pchisq(x, df = df[i], lower.tail = FALSE)
  }
  samp[i] <- x
}

# Computing MLE ----
naive <- mapply(chisq_ncp_mle, x = samp, df = df, normalize = FALSE)
conditional <- mapply(chisq_conditional_ncp_mle, x = samp, df = df,
                      pval_threshold = 0.05, normalize = FALSE)
naivedat <- data.frame(samp = samp, df = df, ncp = ncp,
                       estimate = naive, method = 'naive')
conddat <- data.frame(samp = samp, df = df, ncp = ncp,
                      estimate = conditional, method = "conditional")
plotdat <- rbind(naivedat, conddat)

ggplot(plotdat) +
  geom_point(aes(x = samp, y = estimate^2, col = factor(method))) +
  facet_wrap(~ df, labeller = "label_both") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0) +
  geom_vline(aes(xintercept = qchisq(1 - pval_threshold, df = df)), linetype = 2)

# Computing confidence intervals
naive <- mapply(chisq_conditional_ncp_ci, x = samp, df = df,
                threshold = 0, confidence_level = 0.9)
conditional <- mapply(chisq_conditional_ncp_ci, x = samp, df = df,
                      pval_threshold = 0.05, confidence_level = 0.9)

# plotting
naive_cover <- numeric(n)
conditional_cover <- numeric(n)
for(i in 1:n) {
  naive_cover[i] <- naive[1, i] <= ncp[i] & naive[2, i] >= ncp[i]
  conditional_cover[i] <- conditional[1, i] <= ncp[i] & conditional[2, i] >= ncp[i]
}
mean(naive_cover)
mean(conditional_cover)

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
  facet_grid(df ~ method, labeller = "label_both") +
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






