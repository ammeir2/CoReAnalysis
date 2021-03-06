% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pearson-correlation.R
\name{correlation_conditional_mle}
\alias{correlation_conditional_mle}
\title{Conditional mle for Pearson's correlation}
\usage{
correlation_conditional_mle(corr, df, threshold = NULL,
  pval_threshold = 0.05)
}
\arguments{
\item{corr}{the observed correaltion}

\item{df}{the size of the sample the correaltion was computed from}

\item{threshold}{an alternative way to specify the selection threshold. If specificed,
then the selection rule is assumed to be \code{abs(x) > threshold}. Takes precedence to
\code{pval_threshold}.}

\item{pval_threshold}{the pvalue corresponding to the correlation (under the null)
must be below this value to be observed.}
}
\description{
Conditional mle for Pearson's correlation
}
