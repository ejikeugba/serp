#' Confidence interval for a fitted serp object
#'
#' Provides the confidence interval of estimates for an object of class \code{serp}.
#'
#' @param object An object of class \code{serp}.
#' @param parm unused argument.
#' @param level significance level.
#' @param ... additional arguments.
#' @return A matrix of the the confidence intervals of fitted model.
#' @seealso
#' \code{\link{serp}}, \code{\link{anova.serp}}, \code{\link{vcov.serp}}
#'
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "parallel", link = "logit",
#'            data = wine)
#' confint(m)
#'
#' @export
#'
confint.serp <- function (object, ..., parm, level = 0.95)
{
  ### parm argument is ignored.
  dots <- list(...)
  mod <- as.list(match.call())
  mod[[1]] <- NULL
  mlist <- list(object, ...)
  mc <- unlist(lapply(mlist, class))
  mclass <- all(mc == "serp")
  if (!mclass) stop("input must be an object of class 'serp'")
  if (length(mc) > 1L) stop("one object at a time allowed", call. = FALSE)
  if (!(level > 0 && level < 1 && length(level) == 1 && is.numeric(level)))
    stop("values of 'level' should lie between 0 and 1", call. = FALSE)
  if (!(missing(parm) || is.null(parm)))
    message("argument 'parm' ignored")
  alpha <- 0.5*(1 - level)
  alpha <- c(alpha, 1 - alpha)
  levf <- levformat(alpha, 3)
  qn <- qnorm(alpha)
  est <- coef(object)
  ser <- summary(object)$coefficient[,2]
  cint <- array(NA, dim = c(length(est), 2L),
                dimnames = list(names(est), levf))
  cint[] <- est + ser %o% qn
  return(cint)
}
