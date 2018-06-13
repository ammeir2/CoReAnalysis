#' Converts Pearson's correaltion to z-statstic via Fisher's Transfomration
#'
#' @param corr the observed correaltion
#'
#' @param samp_size the size of the sample the correaltion was computed from
#'
#' @export
corr_to_z <- function(corr, samp_size) {
  fisher_transfrom <- 0.5 * log((1 + corr) / (1 - corr))
  z_stat <- fisher_transfrom * sqrt(samp_size - 3)
  return(z_stat)
}

#' Converts z-statstic to Pearson's correlation via Fisher's Transfomration
#'
#' Should be used to transform back from z-statistic scale to correlation scale
#'
#' @param z_stat the z-statistic
#'
#' @param samp_size the size of the sample the z statistic was computed from
#'
#' @export
z_to_corr <- function(z_stat, samp_size) {
  inverse_transform <- z_stat / sqrt(samp_size - 3)
  corr <- (exp(2 * inverse_transform) - 1) / (1 + exp(2 * inverse_transform))
  return(corr)
}


#' Conditional mle for Pearson's correlation
#'
#' @param corr the observed correaltion
#'
#' @param df the size of the sample the correaltion was computed from
#'
#' @param pval_threshold the pvalue corresponding to the correlation (under the null)
#' must be below this value to be observed.
#'
#' @param threshold an alternative way to specify the selection threshold. If specificed,
#' then the selection rule is assumed to be \code{abs(x) > threshold}. Takes precedence to
#' \code{pval_threshold}.
#'
#' @export
correlation_conditional_mle <- function(corr, df, threshold = NULL,
                                        pval_threshold = 0.05) {
  samp_size <- df
  # Checking threshold
  if(!is.null(threshold)) {
    if(threshold < 0 | threshold >=  1) {
      stop("threshold must be >= 0 and < 1")
    }
  } else {
    threshold <- z_to_corr(qnorm(1 - pval_threshold / 2), samp_size)
  }

  if(abs(corr) < threshold) {
    stop("observed value must exceed the threshold")
  }

  # Computing MLE via Fisher's transform
  z_stat <- corr_to_z(corr, samp_size)
  mle <- t_conditional_ncp_mle(z_stat, threshold = corr_to_z(threshold, samp_size))
  mle <- z_to_corr(mle, samp_size)
  return(mle)
}

#' Correlation Conditional CI
#'
#' @param corr the observed correaltion
#'
#' @param df the size of the sample the correaltion was computed from
#'
#' @param pval_threshold the pvalue corresponding to the correlation (under the null)
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
correlation_conditional_ci <- function(corr, df, threshold = NULL,
                                        pval_threshold = 0.05,
                                        confidence_level = 0.95) {
  samp_size <- df

  # Checking threshold
  if(!is.null(threshold)) {
    if(threshold < 0 | threshold >=  1) {
      stop("threshold must be >= 0 and < 1")
    }
  } else {
    threshold <- z_to_corr(qnorm(1 - pval_threshold / 2), samp_size)
  }

  if(abs(corr) < threshold) {
    stop("observed value must exceed the threshold")
  }

  z_stat <- corr_to_z(corr, samp_size)
  ci <- t_conditional_ncp_ci(z_stat, threshold = corr_to_z(threshold, samp_size),
                             confidence_level = confidence_level)
  ci <- z_to_corr(ci, samp_size)
  return(ci)
}




