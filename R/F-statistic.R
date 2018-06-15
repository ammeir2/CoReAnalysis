#' Compute the mle of the non-centrality parameter of an F test-statistic conditional on selection
#'
#' @param x the observed value
#'
#' @param df1 the numerator degrees of freedom
#'
#' @param df2 the denominator degrees of freedom
#'
#' @param pval_threshold the pvalue corresponding to the chisq-statistic (under the null)
#' must be below this value to be observed.
#'
#' @param threshold an alternative way to specify the selection threshold. If specificed,
#' then the selection rule is assumed to be \code{x > threshold}. Takes precedence to
#' \code{pval_threshold}.
#'
#' @export
f_conditional_ncp_mle <- function(x, df1 = 1, df2 = Inf,
                                  pval_threshold = 0.05,
                                  threshold = NULL) {
  # Checking threshold
  if(!is.null(threshold)) {
    if(threshold < 0) {
      stop("threshold must be non-negative!")
    }
  } else {
    threshold <- qf(1 - pval_threshold, df1 = df1, df2 = df2)
  }

  if(x < threshold) {
    stop("observed value must exceed threshold")
  }

  # Computing the MLE
  mle <- optimize(f = f_conditional_density,
                  df1 = df1, df2 = df2,
                  x = x, threshold = threshold,
                  interval = c(0, x * df1), maximum = TRUE)$maximum

  return(mle)
}

#' The conditional F density
#'
#' @param ncp the non-centrality parameter
#'
#' @param x the observed value
#'
#' @param df the numerator degrees of freedom
#'
#' @param df2 the denominator degress of freedom
f_conditional_density <- function(ncp, x, df1, df2, threshold) {
  dens <- df(x, df1, df2, ncp, TRUE)
  prob <- -pf(threshold, df1, df2, ncp, FALSE, TRUE)
  # print(round(c(x, threshold, ncp, dens + prob), 3))
  return(dens + prob)
}

#' Conditional test-inversion confidence intervals for F ncps
#'
#' @param x the observed value
#'
#' @param df1 the numerator degrees of freedom
#'
#' @param df2 the denominator degrees of freedom
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
f_conditional_ncp_ci <- function(x, df1, df2,
                                 pval_threshold = 0.05,
                                 threshold = NULL,
                                 confidence_level = 0.95) {
  # Checking threshold
  if(!is.null(threshold)) {
    if(threshold < 0) {
      stop("threshold must be non-negative!")
    }
  } else {
    threshold <- qf(1 - pval_threshold, df1 = df1, df2 = df2)
  }

  if(x < threshold) {
    stop("observed value must exceed threshold")
  }

  # Computing CI -----
  lquant <- 1 - (1 - confidence_level) / 2
  uquant <- 1 - lquant
  step_size <- qf(0.55, df1 = df1, df2 = df2) - qf(0.45, df1 = df1, df2 = df2)

  # Finding initial point
  lci <- numerical_invert_ci(cond_f_cdf, lquant, x = x, df1 = df1, df2 = df2,
                             threshold = threshold, step_size = step_size,
                             lbound = 0)
  uci <- numerical_invert_ci(cond_f_cdf, uquant, x = x, df1 = df1, df2 = df2,
                             threshold = threshold, step_size = step_size,
                             lbound = 0)
  ci <- c(lci, uci)
  return(ci)
}


cond_f_cdf <- function(ncp, x, df1, df2, threshold) {
  p_select <- pf(threshold, df1, df2, ncp, lower.tail = FALSE, log.p = FALSE)
  f_cdf <- pf(x, df1, df2, ncp, lower.tail = TRUE, log.p = FALSE)
  return((f_cdf - 1 + p_select) / p_select)
}


