#' Control parameters for serp fit
#'
#' Default control parameters for 'serp' fit. User-supplied control parameters
#' could be specified in the main function.
#'
#' @usage serp.control(
#'              maxits = 5e01,
#'              eps = 1e-07,
#'              maxpen = 1e05,
#'              trace = 0L,
#'              maxAdjIter = 5e0,
#'              max.half.iter = 1e01,
#'              relTol = 1e-03,
#'              nrFold = 5e0,
#'              cv.seed = 1e01,
#'              grid.length = 5e01,
#'              minP = .Machine$double.eps,
#'              ...)
#' @param maxits the maximum number of Newton's iterations. Default to 100.
#' @param eps threshold value during optimization at which the iteration
#'   routine terminates. In other words, when the reported change in the
#'   log-likelihood goes below this threshold, convergence is achieved.
#' @param maxpen the upper end point of the interval from zero to be searched
#'   for a tuning parameter.
#' @param trace prints the Newton's fitting process at each iteration step.If
#'   0 (default) no information is printed, if 1, 2 or 3 different shades of
#'   information are printed.
#' @param maxAdjIter the maximum allowable number of Newton step adjustment to
#'   forestall an early optimization failure. Defaults to 5.
#' @param max.half.iter the maximum number of iteration step-halfings. Defaults
#'   to 10.
#' @param relTol relative convergence tolerance, defaults to 1e-03. checks
#'   relative changes in the parameter estimates between Newton iterations.
#' @param nrFold the number of k-fold cross validation for the CV tuning
#'   method. Default to k = 5.
#' @param cv.seed single numeric value to change the random seed in CV
#'   tuning.
#' @param grid.length the length of the discrete lambda grid for the penalty
#'   method.
#' @param minP A near zero minimum value the fitted probabilities are allowed
#'   to get during iteration to prevent numerical instability .
#' @param ... additional arguments.
#' @return a list of control parameters.
#' @seealso \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#'
#' @export
#'
serp.control <- function(
  maxits = 5e01,
  eps = 1e-07,
  maxpen = 1e05,
  trace = 0L,
  maxAdjIter = 5e0,
  max.half.iter = 1e01,
  relTol = 1e-03,
  nrFold = 5e0,
  cv.seed = 1e01,
  grid.length = 5e01,
  minP = .Machine$double.eps, ...)
{
  cntr <- c(maxits, eps, maxpen, minP, maxAdjIter, trace,
            max.half.iter, relTol)
  if(!all(is.numeric(cntr)) || any((cntr) < 0))
    stop("maxits, eps, maxpen, minP, maxAdjIter, max.half.iter and
  relTol should all be numeric and non-negative")
  if(!(is.numeric(nrFold)) || (nrFold < 2L) || (nrFold > 10L))
    stop("nrFold should be numeric and between 2 and 10 inclusive.",
         call. = FALSE)
  stopcrit <- expression(iter < maxits && sqrt(obs) * abs(objOld - obj)/
                           (abs(objOld) + eps) < eps)
  int.lambdaGrid <- c(10^seq(5, -2, length.out=grid.length), 0)
  msg <- structure(list(
    "0" = "Absolute and relative convergence criteria satisfied",
    "1" = "Absolute convergence criterion satisfied, but relative, criterion
    was not satisfied",
    "2" = "Maximum iteration limit reached",
    "3" = "Maximum number of step-half iteration limit reached",
    "4" = "Maximum number of Newton adjustments reached",
    "5" = "Iteration process abruptly ended",
    "s" = "  Stochastic ordering assumption failed.
    Consider using the penalized, parallel or partial slope,
    or other link functions."
  ))
  list(maxits = as.integer(maxits), eps = eps, minP = minP,
       maxpen = maxpen, maxAdjIter = as.integer(maxAdjIter),
       nrFold = as.integer(nrFold), stopcrit = stopcrit,
       max.half.iter = as.integer(max.half.iter), msg = msg,
       trace = as.integer(trace), cv.seed = cv.seed,
       grid.length = grid.length, relTol = relTol,
       int.lambdaGrid = int.lambdaGrid)
}
