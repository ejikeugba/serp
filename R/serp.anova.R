#' ANOVA method for a fitted serp object
#'
#' Provides a likelihood ratio test for comparing two or more \code{serp} objects. This
#' does not currently support model(s) with penalized slope.
#'
#' @param object An object of class \code{serp}.
#' @param ... additional arguments.
#' @param test type of test to be conducted.
#' @return An ANOVA table with the following components on display:
#' \describe{
#'   \item{model}{the respective model aliases.}
#'   \item{slope}{type of slope fitted, which may be any of, unparallel, parallel,
#'         or partial slope.}
#'   \item{no.par}{the no of parameters in the model.}
#'   \item{AIC}{the akaike information criterion.}
#'   \item{logLik}{the realized log-likelihood.}
#'   \item{Test}{the different pair(s) of test(s) conducted.}
#'   \item{LR.stat}{the computed Likelihood ratio statistic.}
#'   \item{df}{the degree of freedom.}
#'   \item{Pr(chi)}{the p-value of test statitic.}
#' }
#'
#' @seealso
#' \code{\link{serp}}, \code{\link{confint.serp}}, \code{\link{vcov.serp}},
#' \code{\link{errorMetrics}}
#'
#' @examples
#' library(serp)
#' m1 <- serp(rating ~ temp + contact, slope = "parallel", link = "logit",
#'            data = wine)
#' m2 <- update(m1, ~ contact)
#' anova(m1, m2)
#'
#' @export
#'
#'
#'
anova.serp <- function (object, ..., test = c("Chisq", "none"))
{
  mc <- match.call()
  test <- match.arg(test)
  dots <- list(...)
  mod <- as.list(match.call())
  mod[[1]] <- NULL
  nm <- as.character(mod)
  if (!length(unique(nm)) == length(nm))
    stop("duplicate object names are not allowed", call. = FALSE)
  if (!length(dots) && inherits(object, "serp"))
    stop("no anova implementation yet for a single ", "\'serp\' object",
         call. = FALSE)
  mlist <- list(object, ...)
  mclass <- all(unlist(lapply(mlist, class)) == "serp")
  if (!mclass) stop("object(s) of class other than serp not supported",
                    call. = FALSE)
  ml <- length(mlist)
  dr <- vapply(mlist, function(r) (r$nobs - length(r$coef)), numeric(1))
  slope.type <- vapply(mlist, function(r) (r$slope), character(1))
  if (any(slope.type=="penalize"))
    stop("available anova-test doesn't support model(s) with penalized slope",
         call. = FALSE)
  hh <- order(dr, decreasing = TRUE)
  mlist <- mlist[hh]
  nm <- nm[hh]
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
  res <- data.frame(Model = nm, slope = slope.type, no.par =npar, AIC=aic,
                    logLik=llik, Test = prs, LRtest = ch, df = df, Prob = pv)
  names(res) <- c("Model", "slope","no.par", "AIC", "logLik",
                  "Test", "LR.stat", "df", "Pr(Chi)")
  if (test == "none") res <- res[, -9L]
  class(res) <- c("Anova", "data.frame")
  attr(res, "heading") <- "Likelihood ratio tests of ordinal models:\n"
  res
}
