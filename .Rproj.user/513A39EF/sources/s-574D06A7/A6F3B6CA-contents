library(mvtnorm)
library(ggplot2)

sample_corr <- function(n, samp_size, rho) {
  if(length(rho) > 1) {
    stop("only one value of rho can be specified at a time")
  }
  sigma <- matrix(rho, nrow = 2, ncol = 2)
  diag(sigma) <- 1
  samp <- numeric(n)
  norm_samp <- rmvnorm(samp_size * n, sigma = sigma)
  for(i in 1:n) {
    indices <- (samp_size * (i - 1) + 1):(samp_size * i)
    samp[i] <- cor(norm_samp[indices, ])[1, 2]
  }

  return(samp)
}

# Sampling conditionally on selection -----
n <- 1000
samp_size <- sample(c(20, 50, 100), n, replace = TRUE)
rho <- runif(n)
pval_threshold <- 0.05

samp <- numeric(n)
for(i in 1:n) {
  pval <- 1
  while(pval > pval_threshold) {
    corr <- sample_corr(1, samp_size[i], rho[i])
    z <- corr_to_z(corr, samp_size = samp_size[i])
    pval <- 2 * pnorm(-abs(z))
  }
  samp[i] <- corr
}

# Computing MLE
naive <- mapply(correlation_conditional_mle, corr = samp, df = samp_size,
                threshold = 0)
conditional <- mapply(correlation_conditional_mle, corr = samp, df = samp_size,
                      pval_threshold = 0.05)
naivedat <- data.frame(samp = samp, samp_size = samp_size, rho = rho,
                       estimate = naive, method = 'naive')
conddat <- data.frame(samp = samp, samp_size = samp_size, rho = rho,
                      estimate = conditional, method = "conditional")
plotdat <- rbind(naivedat, conddat)

ggplot(plotdat) +
  geom_point(aes(x = samp, y = estimate, col = factor(method))) +
  theme_bw() +
  facet_wrap(~ samp_size, labeller = "label_both") +
  geom_abline(slope = 1, intercept = 0)

# Computing CI
naive <- mapply(correlation_conditional_ci, corr = samp, df = samp_size,
                threshold = 0,
                confidence_level = 0.95)
conditional <- mapply(correlation_conditional_ci, corr = samp, df = samp_size,
                      pval_threshold = .05,
                      confidence_level = 0.95)

naive_cover <- numeric(n)
conditional_cover <- numeric(n)
for(i in 1:n) {
  naive_cover[i] <- naive[1, i] < rho[i] & naive[2, i] > rho[i]
  conditional_cover[i] <- conditional[1, i] < rho[i] & conditional[2, i] > rho[i]
}
mean(naive_cover)
mean(conditional_cover)

naivedat <- data.frame(samp = samp, df = samp_size, rho = rho,
                       cover = naive_cover,
                       lci = naive[1, ],
                       uci = naive[2, ],
                       method = "naive")
conddat <- data.frame(samp = samp, df = samp_size, rho = rho,
                      cover = conditional_cover,
                      lci = conditional[1, ],
                      uci = conditional[2, ],
                      method = "conditional")
plotdat <- rbind(conddat, naivedat)

ggplot(plotdat) +
  geom_segment(aes(y = lci, yend = uci, x = samp, xend = samp, col = factor(cover))) +
  # geom_point(aes(x = samp, y = ncp, fill = df), size = 0.5, shape = 21) +
  theme_bw() +
  ylab("CI") +
  facet_grid(method ~ df, labeller = "label_both") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_abline(slope = 1, intercept = 0, linetype = 2)

ggplot(plotdat) +
  geom_segment(aes(y = lci - rho, yend = uci - rho, x = samp, xend = samp, col = factor(cover))) +
  # geom_point(aes(x = samp, y = ncp, fill = df), size = 0.5, shape = 21) +
  theme_bw() +
  facet_grid(method ~ df, labeller = "label_both") +
  geom_hline(yintercept = 0) +
  ylab("CI - true_ncp")

# Wrapper function -----
obs <- 22
x <- samp[obs]
obs_df <- samp_size[obs]
pval_threshold <- 0.05
fit <- CoReAnalysis(x, test_statistic = "correlation", pval_threshold = 0.05,
                    df1 = obs_df, confidence_level = 1 - 0.5)
replication_power(fit, plot = TRUE)
replication_power(fit, plot = FALSE)






