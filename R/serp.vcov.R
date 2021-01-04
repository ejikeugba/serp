#' Variance covariance matrix for a serp object
#'
#' Provides the Variance covariance matrix of an object of class serp.
#'
#' @param object An object of class serp.
#' @param ... additional arguments.
#' @seealso
#' \code{\link{serp}}
#' @return A variance covariance matrix of a fitted model.
#' @examples
#' \dontrun{
#' vcov(object,...)
#' }
#' @export
#'
vcov.serp <- function(object, ...){
  if (!inherits(object, "serp"))
    stop("not a \"serp\" object", call. = FALSE)
  H <- cbind(object$hess[,c(1:ncol(object$hess))])
  cholHx <- try(chol(H), silent = TRUE)
  if (inherits(cholHx, "try-error"))
    stop("Cannot compute vcov: \nHessian is not positive definite")
  covx   <- chol2inv(cholHx)
  dimnames(covx) <- list(names(object$coef), names(object$coef))
  covx
}
