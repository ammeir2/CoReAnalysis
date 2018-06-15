#' Conditional replication analysis
#'
#' Computes conditional maximum likelihood estimators and confidence
#' confidence intervals for test-statistics that were selected based
#' on their significance.
#'
#' @param x the observed test-statistic.
#'
#' @param test_statistic the type of test-statistic
#'
#' @param pval_threshold the threshold a pvalue must cross (be lesser than)
#' in order for the test-statistic to be observed.
#'
#' @param threshold an alternative way to specifiy a threshold. If a
#' z, t, or correlation then selection rule is assumed to be
#' \code{abs(x) > threshold}. If F or chisq then the selection rule
#' is assumed to be \code{x > threshold}. Takes precedence to
#' \code{pval_threshold}.
#'
#' @param df1 the degrees of freedom for t and chisq statistics. The sample
#' size for Pearson's correlation. The numerator degrees of freedom for an
#' F statistic.
#'
#' @param df2 the denominator degrees of freedom for an F statistic.
#'
#' @param confidence_level the desired confidence level for the
#' confidence interval.
#'
#' @export
CoReAnalysis <- function(x, test_statistic = c("z", "t", "F", "chisq", "correlation"),
                         pval_threshold = 0.05, threshold = NULL,
                         df1 = NULL, df2 = NULL,
                         confidence_level = 0.95, normalize = FALSE) {
  # Checking parameters ----
  check_coreAnalysis_arguments(x, test_statistic, pval_threshold, threshold,
                               df1, df2)
  if(test_statistic == "z") {
    df1 <- Inf
  }

  # Converting thresholds ----
  if(!is.null(threshold)) {
    pval_threshold <- threshold_conversion(x, test_statistic,
                                           threshold = threshold,
                                           pval_threshold = NULL,
                                           df1 = df1, df2 = df2)$pval_threshold
  } else {
    threshold <- threshold_conversion(x, test_statistic,
                                      threshold = NULL,
                                      pval_threshold = pval_threshold,
                                      df1 = df1, df2 = df2)$threshold
  }

  # Computing MLE and CIs-----------
  if(test_statistic %in% c("z", "t")) {
    naive_pvalue <- 2 * pt(-abs(x), df = df1)
    naive_mle <- t_conditional_ncp_mle(x = x, df = df1,
                                       threshold = 0)
    conditional_mle <- t_conditional_ncp_mle(x = x, df = df1,
                                             threshold = threshold)
    naive_ci <- t_conditional_ncp_ci(x, df = df1, threshold = 0,
                                     confidence_level = confidence_level)
    conditional_ci <- t_conditional_ncp_ci(x, df = df1, threshold = threshold,
                                           confidence_level = confidence_level)
  } else if(test_statistic == "chisq") {
    naive_pvalue <- pchisq(x, df = df1, lower.tail = FALSE)
    naive_mle <- chisq_conditional_ncp_mle(x = x, df = df1,
                                       threshold = 0)
    conditional_mle <- t_conditional_ncp_mle(x = x, df = df1,
                                             threshold = threshold)
    naive_ci <- chisq_conditional_ncp_ci(x = x, df = df1, threshold = 0,
                                         confidence_level = confidence_level)
    conditional_ci <- chisq_conditional_ncp_ci(x = x, df = df1,
                                               threshold = threshold,
                                               confidence_level = confidence_level)
  } else if(test_statistic == "correlation") {
    naive_pvalue <- 2 * pnorm(-abs(corr_to_z(x, df1)))
    naive_mle <- correlation_conditional_mle(corr = x, df = df1, threshold = 0)
    conditional_mle <- correlation_conditional_mle(corr = x, df = df1,
                                                   threshold = threshold)
    naive_ci <- correlation_conditional_ci(corr = x, df = df1, threshold = 0,
                                           confidence_level = confidence_level)
    conditional_ci <- correlation_conditional_ci(corr = x, df = df1,
                                                 threshold = threshold,
                                                 confidence_level = confidence_level)
  } else if(test_statistic %in% c("f", "F")) {
    naive_pvalue <- pf(x, df1 = df1, df2 = df2, lower.tail = FALSE)
    naive_mle <- f_conditional_ncp_mle(x, df1 = df1, df2 = df2, threshold = 0)
    conditional_mle <- f_conditional_ncp_mle(x, df1 = df1, df2 = df2,
                                             threshold = threshold)
    conditional_ci <- f_conditional_ncp_ci(x = x, df1 = df1, df2 = df2,
                                           threshold = threshold,
                                           confidence_level = confidence_level)
    naive_ci <- f_conditional_ncp_ci(x = x, df1 = df1, df2 = df2,
                                 threshold = 0,
                                 confidence_level = confidence_level)
  }

  # Normalizing effect size
  if(normalize & test_statistic %in% c("F", "chisq")) {
    naive_mle <- sqrt(naive_mle / df1)
    conditional_mle <- sqrt(conditional_mle / df1)
    naive_ci <- sqrt(naive_ci / df1)
    conditional_ci <- sqrt(conditional_ci / df1)
  }

  # Returning results -------
  result <- list()
  result$x <- x
  result$test_type <- test_statistic
  result$naive_pvalue <- naive_pvalue
  result$threshold <- threshold
  result$pval_threshold <- pval_threshold
  result$df1 <- df1
  result$df2 <- df2
  result$confidence_level <- confidence_level
  result$naive_mle <- naive_mle
  result$conditional_mle <- conditional_mle
  result$naive_ci <- naive_ci
  result$conditional_ci <- conditional_ci
  class(result) <- "CoReAnalysis"
  return(result)
}

#' Pval_threshold to threshold conversion
#'
#' Converts p-value threshold to test-statistic threshold, or
#' the other way around.
#'
#' @param x the observed test-statistic.
#'
#' @param test_statistic the type of test-statistic
#'
#' @param pval_threshold the threshold a pvalue must cross (be lesser than)
#' in order for the test-statistic to be observed. Only one of
#' \code{pval_threshold} or \code{threshold} must be specified, and the other one
#' must be set to NULL.
#'
#' @param threshold an alternative way to specifiy a threshold. If a
#' z, t, or correlation then selection rule is assumed to be
#' \code{abs(x) > threshold}. If F or chisq then the selection rule
#' is assumed to be \code{x > threshold}. Takes precedence to
#' \code{pval_threshold}. Only one of \code{pval_threshold} or
#' \code{threshold} must be specified, and the other one
#' must be set to NULL.
#'
#' @param df1 the degrees of freedom for t and chisq statistics. The sample
#' size for Pearson's correlation. The numerator degrees of freedom for an
#' F statistic.
#'
#' @param df2 the denominator degrees of freedom for an F statistic.
#'
#' @export
threshold_conversion <- function(x, test_statistic = c("z", "t", "F", "chisq", "correlation"),
                                 pval_threshold = 0.05, threshold = NULL,
                                 df1 = NULL, df2 = NULL) {
  # Checking parameters ----
  check_coreAnalysis_arguments(x, test_statistic, pval_threshold, threshold,
                               df1, df2)
  if(test_statistic == "z") {
    df1 <- Inf
  }

  # Checking thresholds ---
  if((is.null(pval_threshold) & is.null(threshold)) |
     (!is.null(pval_threshold) & !is.null(threshold))) {
       stop("exactly one of threshold and pval_threshold must be specified!
            (the other must be null)")
  }

  # Convert threshold to pval_threshold ---
  if(is.null(pval_threshold)) {
    if(test_statistic %in% c("t", "z")) {
      pval_threshold <- pt(-abs(threshold), df = df1) * 2
    } else if(test_statistic == "chisq") {
      pval_threshold <- pchisq(threshold, df = df1, lower.tail = FALSE)
    } else if(test_statistic == "correlation") {
      z <- corr_to_z(threshold, df1)
      pval_threshold <- pnorm(-abs(threshold)) * 2
    } else if(test_statistic %in% c("F", "f")) {
      pval_threshold <- pf(threshold, df1 = df1, df2 = df2, lower.tail = FALSE)
    }
  }

  # Convert pval_threshold to threshold ---
  if(is.null(threshold)) {
    if(test_statistic %in% c("t", "z")) {
      threshold <- qt(1 - pval_threshold / 2, df = df1)
    } else if(test_statistic == "chisq") {
      threshold <- qchisq(1 - pval_threshold, df = df1)
    } else if(test_statistic == "correlation") {
      z_threshold <- qnorm(1 - pval_threshold / 2)
      threshold <- z_to_corr(z_threshold, df1)
    } else if(test_statistic %in% c("F", "f")) {
      threshold <- qf(1 - pval_threshold, df1 = df1, df2 = df2)
    }
  }

  return(list(threshold = threshold, pval_threshold = pval_threshold))
}

#' A helper function for checking coreAnalysis argumnets
check_coreAnalysis_arguments <- function(x, test_statistic = c("z", "t", "F", "chisq", "correlation"),
                                         pval_threshold = 0.05, threshold = NULL,
                                         df1 = NULL, df2 = NULL) {
  if(length(test_statistic) > 1) {
    stop("only a single test_statistic must be specified!")
  }

  if(test_statistic %in% c("t", "chisq", "correlation")) {
    if(is.null(df1)) {
      stop("df1 must specified for t, chisq, and correlation statistics!")
    }
  } else if(test_statistic %in% c("F", "f")) {
    if(is.null(df1) | is.null(df2)) {
      stop("df1 and df2 must specified for F statistics!")
    }

  } else if(test_statistic == "z") {
    df1 <- Inf
  } else {
    stop("specific test statistic not supported!")
  }
}
