#' A function for inverting tests numerically
#'
#' @param cdf a cdf function, the first argument of which is the parameter
#' for which a ci is to be computed.
#'
#' @param quantile the ci quantile to be computed.
#'
#' @param init a value from which to start the search.
#'
#' @param step_size step size for searching for the neighborhood of the correct
#' ci quantile.
#'
#' @param lbound lower bound on the values the parameter can take.
#'
#' @param ubound upper bound on the values the parameter can take.
#'
#' @param ... additional arguments for \code{cdf} function.
numerical_invert_ci <- function(cdf, quantile,
                                init = 0, step_size = 0.1,
                                lbound = -Inf, ubound = Inf,  ...) {
  # bounding the quantile
  lower <- init
  upper <- init
  current <- init
  current_quant <- cdf(current, ...)
  if(current_quant < quantile) {
    while(current_quant < quantile) {
      lower <- lower - step_size
      if(current <= lbound) {
        lower <- lbound
        if(current_quant < quantile) {
          return(lower)
        } else {
          break
        }
      }
      current_quant <- cdf(lower, ...)
      if(is.na(current_quant)) {
        return(lower + step_size)
      }
    }
  } else {
    while(current_quant >= quantile) {
      upper <- upper + step_size
      if(current >= ubound) {
        upper <- ubound
        if(current_quant >= quantile) {
          return(upper)
        } else {
          break
        }
      }
      current_quant <- cdf(upper, ...)
      if(is.na(current_quant)) {
        return(upper - step_size)
      }
    }
  }

  # Searching the quantile on an interval
  ci_quantile <- uniroot(f = function(x) cdf(x, ...) - quantile,
                         interval = c(lower, upper))$root
  return(ci_quantile)
}
