#' Confidence interval for a serp object
#'
#' Provides the confidence interval of estimates of an object of class serp.
#'
#' @param object An object of class serp.
#' @param parm unused argument.
#' @param level significance level.
#' @param ... additional arguments.
#' @return The confidence intervals of a fitted model.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' \dontrun{
#' confint(object,...)
#' }
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
  if (!mclass) stop("input must be an object(s) of class 'serp'")
  if (length(mc) > 1L) stop("one object at a time allowed", call. = FALSE)
  if (!(level > 0 && level < 1 && length(level) == 1 && is.numeric(level)))
    stop("\"clevel!\" should lie between 0 and 1", call. = FALSE)
  if (!(missing(parm) || is.null(parm)))
    message("argument 'parm' ignored")
  alpha <- 0.5*(1 - level)
  alpha <- c(alpha, 1 - alpha)
  levf <- levformat(alpha, 3)
  qn <- qnorm(alpha)
  est <- coef(object)
  ser <- summary(object)$coef[,2]
  cint <- array(NA, dim = c(length(est), 2L),
                dimnames = list(names(est), levf))
  cint[] <- est + ser %o% qn
  return(cint)
}
