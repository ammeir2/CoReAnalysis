#' Power analysis following Conditional Replicability Analysis
#'
#' Computes conditional estimates of power to replicate a study based on
#' a \code{\link[CoReAnalysis]{CoReAnalysis}} fit.
#'
#' @param object a \code{\link[CoReAnalysis]{CoReAnalysis}} object.
#'
#' @param samp_size_ratio the ratio of the sample size in the replication study to
#' the sample size in the original study. For example, if the sample,
#' size in both studies is the same
#' then the ratio is 1 and if the replication study will have 200 subjects where
#' the original study had 100 subjects then the ration will be 100.
#'
#' @param plot if TRUE, then a \code{\link[ggplot2]{ggplot2}} object will be returned.
#'
#' @import ggplot2
#'
#' @export
replication_power <- function(object,
                              samp_size_ratio = 2^(seq(from = -2, to = 3, by = 0.5)),
                              plot = FALSE) {
  if(class(object) != "CoReAnalysis") {
    stop("object must be of class CoReAnalysis")
  }

  x <- object$x
  test_statistic <- object$test_type
  naive_orig_ncp <- object$naive_mle
  conditional_orig_ncp <- object$conditional_mle
  naive_orig_ci <- object$naive_ci
  conditional_orig_ci <- object$conditional_ci
  pval_threshold <- object$pval_threshold

  # Computing power estimates ----
  if(test_statistic %in% c("t", "z")) {
    naive_orig_ci <- sort(naive_orig_ci * sign(x))
    if(sign(naive_orig_ci[1]) != sign(naive_orig_ci[2])) {
      naive_orig_ci[1] <- 0
    }
    conditional_orig_ci <- sort(conditional_orig_ci * sign(x))
    if(sign(conditional_orig_ci[1]) != sign(conditional_orig_ci[2])) {
      conditional_orig_ci[1] <- 0
    }

    # ncp point estimates
    df <- object$df1
    rep_df <- df * samp_size_ratio
    naive_rep_ncp <- naive_orig_ncp * sqrt(samp_size_ratio)
    conditional_rep_ncp <- conditional_orig_ncp * sqrt(samp_size_ratio)

    # ncp cis
    naive_rep_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    conditional_rep_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    naive_rep_ci[, 1] <- naive_orig_ci[1] * sqrt(samp_size_ratio)
    naive_rep_ci[, 2] <- naive_orig_ci[2] * sqrt(samp_size_ratio)
    conditional_rep_ci[, 1] <- conditional_orig_ci[1] * sqrt(samp_size_ratio)
    conditional_rep_ci[, 2] <- conditional_orig_ci[2] * sqrt(samp_size_ratio)

    # power point estimates
    rep_threshold <- qt(1 - pval_threshold / 2, df = rep_df)
    naive_power <- two_sided_t(rep_threshold, rep_df, naive_rep_ncp)
    conditional_power <- two_sided_t(rep_threshold, rep_df, conditional_rep_ncp)

    # power CIs
    naive_power_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    naive_power_ci[, 1] <- two_sided_t(rep_threshold, rep_df, naive_rep_ci[, 1])
    naive_power_ci[, 2] <- two_sided_t(rep_threshold, rep_df, naive_rep_ci[, 2])
    conditional_power_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    conditional_power_ci[, 1] <- two_sided_t(rep_threshold, rep_df, conditional_rep_ci[, 1])
    conditional_power_ci[, 2] <- two_sided_t(rep_threshold, rep_df, conditional_rep_ci[, 2])

    # ncp CIs with no sign correction
    naive_orig_ci <- object$naive_ci
    conditional_orig_ci <- object$conditional_ci
    naive_rep_ci[, 1] <- naive_orig_ci[1] * sqrt(samp_size_ratio)
    naive_rep_ci[, 2] <- naive_orig_ci[2] * sqrt(samp_size_ratio)
    conditional_rep_ci[, 1] <- conditional_orig_ci[1] * sqrt(samp_size_ratio)
    conditional_rep_ci[, 2] <- conditional_orig_ci[2] * sqrt(samp_size_ratio)
  }
  else if(test_statistic == "correlation") {
    naive_orig_ci <- sort(naive_orig_ci * sign(x))
    if(sign(naive_orig_ci[1]) != sign(naive_orig_ci[2])) {
      naive_orig_ci[1] <- 0
    }
    conditional_orig_ci <- sort(conditional_orig_ci * sign(x))
    if(sign(conditional_orig_ci[1]) != sign(conditional_orig_ci[2])) {
      conditional_orig_ci[1] <- 0
    }

    # point estimates
    naive_rho <- object$naive_mle
    cond_rho <- object$conditional_mle
    samp_size <- object$df1

    # ncp CIs don't change with the sample size for correlations
    # Converting to z-values for simpler calculations
    naive_ncp <- corr_to_z(naive_rho, samp_size) * sqrt(samp_size_ratio)
    cond_ncp <- corr_to_z(cond_rho, samp_size) * sqrt(samp_size_ratio)
    naive_ncp_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    cond_ncp_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    naive_ncp_ci[, 1] <- corr_to_z(naive_orig_ci[1], samp_size) * sqrt(samp_size_ratio)
    naive_ncp_ci[, 2] <- corr_to_z(naive_orig_ci[2], samp_size) * sqrt(samp_size_ratio)
    cond_ncp_ci[, 1] <- corr_to_z(conditional_orig_ci[1], samp_size) * sqrt(samp_size_ratio)
    cond_ncp_ci[, 2] <- corr_to_z(conditional_orig_ci[2], samp_size) * sqrt(samp_size_ratio)
    threshold <- qnorm(1 - pval_threshold / 2)

    # power point estimates
    naive_power <- two_sided_t(threshold, Inf, naive_ncp)
    conditional_power <- two_sided_t(threshold, Inf, cond_ncp)

    # power cis
    naive_power_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    naive_power_ci[, 1] <- two_sided_t(threshold, Inf, naive_ncp_ci[, 1])
    naive_power_ci[, 2] <- two_sided_t(threshold, Inf, naive_ncp_ci[, 2])
    conditional_power_ci <- matrix(nrow = length(samp_size_ratio), ncol = 2)
    conditional_power_ci[, 1] <- two_sided_t(threshold, Inf, cond_ncp_ci[, 1])
    conditional_power_ci[, 2] <- two_sided_t(threshold, Inf, cond_ncp_ci[, 2])

    # ncp CIs with no sign correction
    naive_rep_ci <- object$naive_ci
    conditional_rep_ci <- object$conditional_ci
    naive_rep_ncp <- object$naive_mle
    conditional_rep_ncp <- object$conditional_mle
  }
  else if(test_statistic %in% c("chisq", "f", "F")) {
    df1 <- object$df1
    df2 <- object$df2
    threshold <- object$threshold
    if(test_statistic == "chisq") {
      df2 <- Inf
      threshold <- object$threshold / df1
    }

    # ncp power estimates
    naive_rep_ncp <- naive_orig_ncp * samp_size_ratio
    conditional_rep_ncp <- conditional_orig_ncp * samp_size_ratio

    # ncp CIs
    naive_rep_ci <- matrix(ncol = 2, nrow = length(samp_size_ratio))
    naive_rep_ci[, 1] <- naive_orig_ci[1] * samp_size_ratio
    naive_rep_ci[, 2] <- naive_orig_ci[2] * samp_size_ratio
    conditional_rep_ci <- matrix(ncol = 2, nrow = length(samp_size_ratio))
    conditional_rep_ci[, 1] <- conditional_orig_ci[1] * samp_size_ratio
    conditional_rep_ci[, 2] <- conditional_orig_ci[2] * samp_size_ratio

    # power point estimates
    naive_power <- 1 - pf(threshold, df1, df2, ncp = naive_rep_ncp)
    conditional_power <- 1 - pf(threshold, df1, df2, ncp = conditional_rep_ncp)

    # power CIs
    naive_power_ci <- matrix(ncol = 2, nrow = length(samp_size_ratio))
    naive_power_ci[, 1] <- 1 - pf(threshold, df1, df2, ncp = naive_rep_ci[, 1])
    naive_power_ci[, 2] <- 1 - pf(threshold, df1, df2, ncp = naive_rep_ci[, 2])
    conditional_power_ci <- matrix(ncol = 2, nrow = length(samp_size_ratio))
    conditional_power_ci[, 1] <- 1 - pf(threshold, df1, df2, ncp = conditional_rep_ci[, 1])
    conditional_power_ci[, 2] <- 1 - pf(threshold, df1, df2, ncp = conditional_rep_ci[, 2])
  }

  # plotting -----
  if(plot) {
    naivedat <- data.frame(samp_size_ratio = samp_size_ratio,
                           lci = naive_power_ci[, 1],
                           uci = naive_power_ci[, 2],
                           estimate = naive_power,
                           offset = 0,
                           method = "naive")
    conddat <- data.frame(samp_size_ratio = samp_size_ratio,
                          lci = conditional_power_ci[, 1],
                          uci = conditional_power_ci[, 2],
                          estimate = conditional_power,
                          offset = 0.05,
                          method = "conditional")
    plotdat <- rbind(naivedat, conddat)
    plotdat <- plotdat[order(plotdat$samp_size_ratio, decreasing = FALSE), ]
    plot_object <- ggplot(plotdat, aes(x = log2(samp_size_ratio) + offset, col = method, linetype = method, shape = method)) +
      geom_point(aes(y = estimate)) +
      geom_segment(aes(xend = log2(samp_size_ratio) + offset, y = lci, yend = uci)) +
      theme_bw() +
      # facet_wrap(~ method, labeller = "label_both") +
      ylab("Power Estimate") +
      xlab("log2(replication sample size/original sample size)")
    return(plot_object)
  }

  # Reporting -----
  object$samp_size_ratio <- samp_size_ratio

  # ncp
  object$naive_rep_ncp <- naive_rep_ncp
  object$conditional_rep_ncp <- conditional_rep_ncp
  object$naive_rep_ncp_ci <- naive_rep_ci
  object$conditional_rep_ncp_ci <- conditional_rep_ci

  # power
  object$naive_power <- naive_power
  object$conditional_power <- conditional_power
  object$naive_power_ci <- naive_power_ci
  object$conditional_power_ci <- conditional_power_ci

  # fin
  class(object) <- "CoRe_power"
  return(object)
}

#' A helper function for computing the selection probability of a t test-statistic
#'
#' @param threshold the threhold that must be crossed (in absolute value).
#'
#' @param df the degrees of freedom.
#'
#' @param the non-centrality parameter.
#'
two_sided_t <- function(threshold, df, ncp) {
  1 + pt(-threshold, df, ncp) - pt(threshold, df, ncp)
}
