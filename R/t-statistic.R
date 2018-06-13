#' A function for estimating the MLE of the ncp of a t distribution
#'
#' @param x the observed value
#'
#' @param df the degrees of freedom of \code{x}.
#'
#' @export
t_ncp_mle <- function(x, df = Inf) {
  mle <- nlm(f = function(ncp) -dt(x, df = df, ncp, log = TRUE), p = x)$estimate
  return(mle)
}

#' Computes a naive CI for the t non-centrality parameter
#'
#' @param x the observed value
#'
#' @param df the degrees of freedom of \code{x}.
#'
#' @param confidence_level the desired confidence level of the confidence
#' interval.
#'
#' @export
t_ncp_ci <- function(x, df = Inf, confidence_level = 0.95) {
  lquant <- 1 - (1 - confidence_level) / 2
  uquant <- 1 - lquant
  tcdf <- function(ncp, x, df) pt(x, df, ncp)
  lci <- numerical_invert_ci(cdf, lquant, x = x, df = df)
  uci <- numerical_invert_ci(cdf, uquant, x = x, df = df)
  return(c(lci, uci))
}

#' A function for computing the MLE of a t distribution conditional on selection
#'
#' @param x the observed value
#'
#' @param df the degrees of freedom of \code{x}.
#'
#' @param pval_threshold the pvalue corresponding to the t-statistic (under the null)
#' must be below this value to be observed.
#'
#' @param threshold an alternative way to specify the selection threshold. If specificed,
#' then the selection rule is assumed to be \code{abs(x) > threshold}. Takes precedence to
#' \code{pval_threshold}.
#'
#' @param iterlim maximum number of numerical optimization iterations.
#' See \code{\link[stats]{nlm}} for details.
#'
#' @export
t_conditional_ncp_mle <- function(x, df = Inf, pval_threshold, threshold = NULL,
                                  iterlim = 1000) {
  # Check the threshold
  if(!is.null(threshold)) {
    if(threshold < 0) {
      stop("threshold must be non-negative!")
    }
  } else {
    threshold <- qt(1 - pval_threshold / 2, df = df)
  }

  # Check the observed value
  if(abs(x) < threshold) {
    stop("Observed value should exceed the threshold!")
  }

  # Computing the MLE
  mle <- nlm(f = t_conditional_density, p = x,
             x = x, df = df, threshold = threshold)$estimate
  return(mle)
}

#' A function for computing the (negative) conditional t density
#'
#' @param x the observed value
#'
#' @param df the degrees of freedom of \code{x}.
#'
#' @param threshold The threshold defines the selection rule, \code{abs(x) > threshold}.
t_conditional_density <- function(ncp, x, df = Inf, threshold) {
  logdens <- dt(x, df, ncp, log = TRUE)
  uprob <- pt(threshold, df = df, ncp = ncp, lower.tail = FALSE, log.p = TRUE)
  lprob <- pt(-threshold, df = df, ncp = ncp, lower.tail = TRUE, log.p = TRUE)
  if(ncp >= 0) {
    logprob <- uprob + log(1 + exp(lprob - uprob))
  } else {
    logprob <- lprob + log(1 + exp(uprob - lprob))
  }
  return(-logdens + logprob)
}

#' Compute a conditional ci for a t-statistic non-centraility parameter
#'
#' @param x the observed value
#'
#' @param df the degrees of freedom of \code{x}.
#'
#' @param pval_threshold the pvalue corresponding to the t-statistic (under the null)
#' must be below this value to be observed.
#'
#' @param threshold an alternative way to specify the selection threshold. If specificed,
#' then the selection rule is assumed to be \code{abs(x) > threshold}. Takes precedence to
#' \code{pval_threshold}.
#'
#' @param confidence_level the desired confidence level of the confidence
#' interval.
#'
#' @export
t_conditional_ncp_ci <- function(x, df = Inf, pval_threshold, threshold = NULL,
                                 confidence_level = 0.95) {
  # Check the threshold
  if(!is.null(threshold)) {
    if(threshold < 0) {
      stop("threshold must be non-negative!")
    }
  } else {
    threshold <- qt(1 - pval_threshold / 2, df = df)
  }

  # Check the observed value
  if(abs(x) < threshold) {
    stop("Observed value should exceed the threshold!")
  }

  # Converting observed value to negative
  if(x >= threshold) {
    x <- -x
    end_sign <- -1
  } else {
    end_sign <- 1
  }

  # Computing CI
  cond_t_cdf <- function(ncp, x, df, threshold) {
    numerator <- pt(x, df, ncp, log.p = TRUE)
    lprob <- pt(-threshold, df, ncp = ncp, log.p = TRUE)
    uprob <- pt(threshold, df, ncp = ncp, log.p = TRUE, lower.tail = FALSE)
    denom <- lprob + log(1 + exp(uprob - lprob))
    return(exp(numerator - denom))
  }

  lquant <- 1 - (1 - confidence_level) / 2
  uquant <- 1 - lquant
  lci <- numerical_invert_ci(cond_t_cdf, lquant, x = x, df = df, threshold = threshold)
  uci <- numerical_invert_ci(cond_t_cdf, uquant, x = x, df = df, threshold = threshold)
  ci <- sort(end_sign * c(lci, uci))
  return(ci)
}



