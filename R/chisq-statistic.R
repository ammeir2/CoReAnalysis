#' Compute the mle of the non-centrality parameter for a chi-square test-statistic
#'
#' @param x the observed value
#'
#' @param df the degrees of freedom
#'
#' @param normalize whether to report the non-centrality parameter per degree of
#' freedom
#'
#' @export
chisq_ncp_mle <- function(x, df, normalize = TRUE) {
  mle <- optimize(f = function(ncp) dchisq(x, df, ncp, log = TRUE),
                  interval = c(0, x), maximum = TRUE)$maximum
  if(normalize) {
    mle <- sqrt(mle) / df
  } else {
    mle <- sqrt(mle)
  }
  return(mle)
}

#' Compute the mle of the non-centrality parameter of a chi-square test-statistic conditional on selection
#'
#' @param x the observed value
#'
#' @param df the degrees of freedom
#'
#' @param pval_threshold the pvalue corresponding to the chisq-statistic (under the null)
#' must be below this value to be observed.
#'
#' @param threshold an alternative way to specify the selection threshold. If specificed,
#' then the selection rule is assumed to be \code{x > threshold}. Takes precedence to
#' \code{pval_threshold}.
#'
#' #' @param normalize whether to report the non-centrality parameter per degree of
#' freedom
#'
#' @export
chisq_conditional_ncp_mle <- function(x, df, pval_threshold = 0.05,
                                      threshold = NULL,
                                      normalize = FALSE) {
  # Checking threshold
  if(!is.null(threshold)) {
    if(threshold < 0) {
      stop("threshold must be non-negative!")
    }
  } else {
    threshold <- qchisq(1 - pval_threshold, df = df)
  }

  if(x < threshold) {
    stop("observed value must exceed threshold")
  }

  # Computing the MLE
  mle <- optimize(f = chisq_conditional_density,
                  df = df, x = x, threshold = threshold,
                  interval = c(0, x), maximum = TRUE)$maximum

  if(normalize) {
    mle <- sqrt(mle) / df
  } else {
    mle <- sqrt(mle)
  }

  return(mle)
}

#' The conditional chi-square density
#'
#' @param ncp the non-centrality parameter
#'
#' @param x the observed value
#'
#' @param df the chi-square degrees of freedom
chisq_conditional_density <- function(ncp, x, df, threshold) {
  dchisq(x, df, ncp, TRUE) - pchisq(threshold, df, ncp, FALSE, TRUE)
}

#' Conditional test-inversion confidence intervals for chi-square ncps
#' @param x the observed value
#'
#' @param df the degrees of freedom
#'
#' @param pval_threshold the pvalue corresponding to the chisq-statistic (under the null)
#' must be below this value to be observed.
#'
#' @param threshold an alternative way to specify the selection threshold. If specificed,
#' then the selection rule is assumed to be \code{x > threshold}. Takes precedence to
#' \code{pval_threshold}.
#'
#' @param confidence_level the desired confidence level of the confidence
#' interval.
#'
#' @param normalize whether to report the non-centrality parameter per degree of
#' freedom
#'
#' @export
chisq_conditional_ncp_ci <- function(x, df, pval_threshold = 0.05,
                                     threshold = NULL,
                                     confidence_level = 0.95,
                                     normalize = FALSE) {
  # Checking threshold
  if(!is.null(threshold)) {
    if(threshold < 0) {
      stop("threshold must be non-negative!")
    }
  } else {
    threshold <- qchisq(1 - pval_threshold, df = df)
  }

  if(x < threshold) {
    stop("observed value must exceed threshold")
  }

  # Computing CI -----
  cond_chisq_cdf <- function(ncp, x, df, threshold) {
    p_select <- pchisq(threshold, df, ncp, lower.tail = FALSE, log.p = FALSE)
    chisq_cdf <- pchisq(x, df, ncp, lower.tail = TRUE, log.p = FALSE)
    return((chisq_cdf - 1 + p_select) / p_select)
  }

  lquant <- 1 - (1 - confidence_level) / 2
  uquant <- 1 - lquant
  lci <- numerical_invert_ci(cond_chisq_cdf, lquant, x = x, df = df, threshold = threshold,
                             step_size = df, lbound = 0)
  uci <- numerical_invert_ci(cond_chisq_cdf, uquant, x = x, df = df, threshold = threshold,
                             step_size = df, lbound = 0)
  ci <- c(lci, uci)
  return(ci)
}




