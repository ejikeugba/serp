#' ANOVA method for an object of class 'serp'
#'
#' Provides an ANOVA table for comparing two or more 'serp' objects.
#'
#' @param object An object of class 'serp'.
#' @param ... additional arguments.
#' @param test type of test to be conducted.
#' @return An ANOVA table with the following components:
#' \describe{
#'   \item{model}{the respective model aliases.}
#'   \item{no.par}{the no of parameters in the model.}
#'   \item{AIC}{the akaike information criterion.}
#'   \item{logLik}{the realized log-likelihood.}
#'   \item{Test}{the different pair(s) of test(s) conducted.}
#'   \item{LR.stat}{the computed Likelihood ratio statistic.}
#'   \item{df}{the degree of freedom.}
#'   \item{Pr(chi)}{the p-value of test statitic.}
#' }
#'
#' @seealso \code{\link{serp}}
#'
#' @examples
#' # See serp() documentation for examples.
#'
#' @export
#'
anova.serp <- function (object, ..., test = c("Chisq", "none"))
{
  test <- match.arg(test)
  dots <- list(...)
  mod <- as.list(match.call())
  mod[[1]] <- NULL
  nm <- as.character(mod)
  if (!length(dots) && inherits(object, "serp"))
    stop("no anova implementation yet for a single ", "\'serp\' object",
         call. = FALSE)
  mlist <- list(object, ...)
  mclass <- all(unlist(lapply(mlist, class)) == "serp")
  if (!mclass) stop("object(s) not of class \'serp\'")
  ml <- length(mlist)
  dr <- vapply(mlist, function(r) (r$nobs - length(r$coef)), numeric(1))
  hh <- order(dr, decreasing = TRUE)
  mlist <- mlist[hh]
  if (any(!vapply(mlist, inherits, logical(1), "serp")))
    stop("input must be an object(s) of class 'serp'", call. = FALSE)
  fv <- vapply(mlist, function(r) length(r$fitted.values), numeric(1))
  if (any(fv != fv[1L]))
    stop("models were not all fitted to the same ",
         "size of dataset", call. = FALSE)
  dep <- unique(vapply(mlist, function(r) paste(formula(r)[2L]), character(1)))
  drs <- dr[hh]
  dev <- vapply(mlist, function(r) -2*r$logLik, numeric(1))
  aic <- round(vapply(mlist, function(r) r$aic, numeric(1)), 2L)
  npar <- vapply(mlist, function(r) length(r$coef), numeric(1))
  llik <- round(vapply(mlist, function(r) r$logLik, numeric(1)), 3L)
  prs <- c("", paste(seq_len(ml - 1L), 2L:ml, sep = " vs "))
  df <- c(NA_integer_, -diff(drs))
  ch <- round(c(NA_real_, -diff(dev)), 3L)
  pv <- c(NA_real_, 1 - pchisq(ch[-1L], df[-1L]))
  pv <- signif(pv, 4)
  res <- data.frame(Model = nm, no.par =npar, AIC=aic, logLik=llik,
                    Test = prs, LRtest = ch, df = df, Prob = pv)
  names(res) <- c("Model", "no.par", "AIC", "logLik",
                  "Test", "LR.stat", "df", "Pr(Chi)")
  if (test == "none") res <- res[, -7L]
  class(res) <- c("Anova", "data.frame")
  attr(res, "heading") <- c("Likelihood ratio tests of ordinal models.",
                            paste("Response:",dep, "\n"))
  res
}
