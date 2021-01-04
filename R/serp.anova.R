#' An ANOVA method for an object of class serp
#'
#' Provides an ANOVA table for comparing two or more serp objects.
#'
#' @param object An object of class serp.
#' @param ... additional arguments.
#' @param test type of test to be conducted.
#' @return The ANOVA table of a fitted model.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' \dontrun{
#' anova(object,..., test)
#' }
#' @export
#'
anova.serp <- function (object, ..., test = c("Chisq", "none"))
{
  if (!inherits(object, "serp"))
    stop("not a \"serp\" object", call. = FALSE)
  test <- match.arg(test)
  dots <- list(...)
  mod <- as.list(match.call())
  mod[[1]] <- NULL
  nm <- as.character(mod)
  if (!length(dots))
    stop("no anova implementation for a single ",
         "\"serp\" object", call. = FALSE)
  mlist <- list(object, ...)
  ml <- length(mlist)
  dr <- sapply(mlist, function(r) (r$nobs - length(r$coef)))
  hh <- order(dr, decreasing = TRUE)
  mlist <- mlist[hh]
  if (any(!sapply(mlist, inherits, "serp")))
    stop("objects not of class \"serp\"", call. = FALSE)
  fv <- sapply(mlist, function(r) length(r$fitted.values))
  if (any(fv != fv[1L]))
    stop("models were not all fitted to the same ",
         "size of dataset", call. = FALSE)
  dep <- unique(sapply(mlist, function(r) paste(formula(r)[2L])))
  drs <- dr[hh]
  dev <- sapply(mlist, function(r) -2*r$logLik)
  aic <- round(sapply(mlist, function(r) r$aic), 2L)
  npar <- sapply(mlist, function(r) length(r$coef))
  llik <- round(sapply(mlist, function(r) r$logLik), 3L)
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
