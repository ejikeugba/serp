#' Control parameters for serp fit.
#'
#' Default control parameters for serp fit. User supplied control parameters could be
#' specified in the main function.
#'
#' @usage serp.control(maxits = 5e01, eps = 1e-07, maxpen = 1e5,
#'                     minP = .Machine$double.eps, nFold = 5e0)
#' @param maxits the maximum number of iteration.
#' @param eps threshold value during optimization at which the iteration routine terminates. In other words, when the reported change in the log-likelihood goes below this threshold, convergence is achieved.
#' @param maxpen the upper end point of the interval from zero to be searched for a minimum tuning parameter.
#' @param minP A near zero minimum value the fitted probabilities are allowed to get during iteration to prevent numerical instability .
#' @param nFold the number of k-fold cross validation for the "cv" \code{tuning} method, with default value k = 5.
#' @return NULL
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#'
#' @export
#'
serp.control <- function(maxits = 5e01, eps = 1e-07, maxpen = 1e05,
                         minP = .Machine$double.eps, nFold = 5e0)
{
  if (!is.numeric(eps) || eps <= 0)
    stop("value of eps must be > 0")
  if (!is.numeric(maxits) || maxits <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(maxpen) || maxpen <= 0)
    stop("value of maxpen must be > 0")
  if (!is.numeric(minP) || minP <= 0)
    stop("value of minP must be > 0")
  if (is.null(maxits)) maxits <- 5e01
  if (is.null(eps)) eps <- 1e-07
  if (is.null(minP)) minP <- .Machine$double.eps
  if (is.null(maxpen)) maxpen <- 1e05
  if(!(is.numeric(nFold)) || (nFold < 2L) || (nFold > 10L) ||
     (round(nFold) != nFold))
    stop("nFold should be numeric, non-decimal, ",
         "and between 2 and 10 inclusive.", call. = FALSE)
  if (is.null(nFold)) nFold <- 5e0
  stopcrit <- expression(
    iter < maxits && sqrt(obs) * abs(objOld - objNew)/
      (abs(objOld) + eps) < eps)
  list(maxits = maxits, eps = eps, maxpen = maxpen,
       minP=minP, nFold = nFold, stopcrit=stopcrit)
}
