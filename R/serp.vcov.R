#' Variance covariance matrix for a fitted serp object
#'
#' Provides the Variance covariance matrix of an object of class \code{serp}.
#'
#' @param object An object of class \code{serp}.
#' @param ... additional arguments.
#' @seealso \code{\link{serp}}
#' @return A variance covariance matrix of a fitted model.
#' @seealso
#' \code{\link{serp}}, \code{\link{anova.serp}}, \code{\link{confint.serp}}
#'
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "parallel", link = "logit",
#'            data = serp::wine)
#' vcov(m)
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
  if (!mclass) stop("input must be an object of class 'serp'")
  if (length(mc) > 1L) stop("one object at a time allowed", call. = FALSE)
  H <- cbind(object$hess[,seq_len(ncol(object$hess))])
  cholHx <- try(chol(H), silent = TRUE)
  if (inherits(cholHx, "try-error"))
    stop("Could not compute vcov: \nHessian is not positive definite")
  covx   <- chol2inv(cholHx)
  dimnames(covx) <- list(names(object$coef), names(object$coef))
  covx
}
