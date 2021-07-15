#' Variance covariance matrix for a serp object
#'
#' Provides the Variance covariance matrix of an object of class 'serp'.
#'
#' @param object An object of class 'serp'.
#' @param ... additional arguments.
#' @seealso \code{\link{serp}}
#' @return A variance covariance matrix of a fitted model.
#'
#' @examples
#' # See serp() documentation for examples.
#'
#' @export
#'
vcov.serp <- function(object, ...){
  dots <- list(...)
  mod <- as.list(match.call())
  mod[[1]] <- NULL
  mlist <- list(object, ...)
  mc <- unlist(lapply(mlist, class))
  mclass <- all(mc == "serp")
  if (!mclass) stop("input must be an object(s) of class 'serp'")
  if (length(mc) > 1L) stop("one object at a time allowed", call. = FALSE)

  if (!inherits(object, "serp"))
    stop("not a \"serp\" object", call. = FALSE)
  H <- cbind(object$hess[,seq_len(ncol(object$hess))])
  cholHx <- try(chol(H), silent = TRUE)
  if (inherits(cholHx, "try-error"))
    stop("Could not compute vcov: \nHessian is not positive definite")
  covx   <- chol2inv(cholHx)
  dimnames(covx) <- list(names(object$coef), names(object$coef))
  covx
}
