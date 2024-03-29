#' Print method for a fitted serp object
#'
#' Prints out a vector of coefficients of the fitted model with some
#' additional goodness-of-fit measures.
#'
#' @param x An object of class \code{serp}.
#' @param ... additional arguments.
#' @return No return value
#' @seealso \code{\link{serp}}, \code{\link{print.summary.serp}}
#'
#' @export
#'
print.serp <- function(x, ...)
{
  if (!inherits(x, "serp"))
    stop("input must be an object of class 'serp'", call. = FALSE)
  object <- x
  cat(crayon::green("\ncall:\n"))
  print(object$call)
  max.tun <- FALSE
  cat(crayon::green("\ncoefficients:\n"))
  print(object$coef)
  cat("\nloglik:", object$logLik, " ","aic:", object$aic, "\n")

  if (object$slope == 'penalize') max.tun <- penalty.print(object, max.tun)

  if (max.tun){
    cat("\n")
    message("* optimal tuning obtained at lambda grid upper limit\n")
  }
  if (!length(colnames(object$model)[-1L])==0L){
    if (length(colnames(object$model)[-1L]) ==
        length(object$globalEff)) message("\nno subject-specific effects found, ",
                                          "parallel slope returned!\n")
  }
  na.ac <- length(object$na.action)
  if (na.ac > 0){
    cat("\n")
    message("\n--- ",na.ac," observation(s) deleted due to ",
            "missingness ","---", "\n")
  }
  invisible(object)
}

#' Print method for an object of class summary.serp
#'
#' Prints out the information supplied via \code{summary.serp} method.
#'
#' @param x An object of class \code{summary.serp}.
#' @param ... additional arguments.
#' @return No return value
#' @seealso \code{\link{serp}}, \code{\link{print.serp}}
#'
#' @export
#'
print.summary.serp <- function(x, ...){
  if (!inherits(x, "summary.serp"))
    stop("input must be an object of class 'serp'", call. = FALSE)
  object <- x
  cat(crayon::green("\ncall:\n"))
  print(object$call)
  max.tun <- FALSE
  nL <- object$ylev
  coef <- as.data.frame(object$coefficients)
  df <- object$rdf
  na.ac <- length(object$na.action)
  cat(crayon::green("\nCoefficients:\n"))
  printCoefmat(coef, digits = 4L, signif.stars = TRUE,
               na.print = "NA", P.values = TRUE, has.Pvalue = TRUE,
               signif.legend = TRUE)
  cat("\nNumber of iterations:", object$iter, "\n")
  cat("\nLoglik:", object$logLik, "on", df, "degrees of freedom","\n")
  cat("\nAIC:", object$aic)
  cat("\n")
  if (!ncol(object$model) == 1L){
    cat(crayon::green("\nExponentiated coefficients:\n"))
    print(object$expcoefs)
  }
  if (object$slope == 'penalize') max.tun <- penalty.print(object, max.tun)
  if (max.tun){
    cat("\n")
    message("* optimal tuning obtained at lambda grid upper limit\n")
  }
  if (!length(colnames(object$model)[-1L])==0L){
    if (length(colnames(object$model)[-1L]) ==
        length(object$globalEff)) message("\nno subject-specific effects found, ",
                                          "parallel slope returned!\n")
  }
    if (na.ac > 0){
    cat("\n")
    message("\n--- ",na.ac," observation(s) deleted due to ",
            "missingness ","---", "\n")
  }
  invisible(object)
}

#' Summary method for a fitted serp object.
#'
#' This function summarizes the result of a fitted serp object in a dataframe.
#'
#' @param object An object of class \code{serp}.
#' @param ... Not used. Additional summary arguments.
#' @details
#' an object of class \code{summary.serp}. A list (depending on the type of
#' \code{slope} used) of all model components defined in the \code{\link{serp}},
#' function with additional components listed below.
#' @return \item{coefficients}{the matrix of coefficients, standard errors,
#'         z-values and p-values.}
#' @return \item{null.deviance}{the deviance for the intercept only model.}
#' @return \item{null.logLik}{the log-likelihood for the intercept only model.}
#' @return \item{penalty}{list of penalization information obtained with
#'         \code{slope} set to "penalize".}
#' @return \item{expcoefs}{the exponentiated coefficients.}
#'
#' @seealso \code{\link{anova.serp}}, \code{\link{predict.serp}},
#' \code{\link{confint.serp}}, \code{\link{vcov.serp}}
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "penalize",
#'            reverse = TRUE, link = "logit", tuneMethod = "user",
#'            lambda = 0, data = wine)
#' summary(m)
#' @export
#'
summary.serp <- function(object, ...){
  if (!inherits(object, "serp"))
    stop("input must be an object of class 'serp'", call. = FALSE)
  coef <- object$coef
  nL <- object$ylev
  coefs <- matrix(0, length(coef), 4L,
                  dimnames = list(names(coef), NULL))
  coefs[, 1L] <- coef
  H <- cbind(object$hessian[,seq_len(ncol(object$hessian))])
  dimnames(H) <- list(names(object$coef), names(object$coef))
  cholHx <- try(chol(H), silent = TRUE)
  if (inherits(cholHx, "try-error"))
    cholHx <- covx <- diag(NA, dim(H))
  else covx   <- chol2inv(cholHx)
  se.est <- sqrt(diag(covx))
  coefs[,2L] <- se.est
  coefs[,3L] <- z.value <- coef/se.est
  coefs[,4L] <- pvalue  <- 2 * pnorm(abs(z.value), lower.tail = FALSE)
  colnames(coefs) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)")
  expcoefs <- exp(coef(object)[-c(1:(nL-1))])
  if (object$slope == 'penalize'){
    tun <- object$tuneMethod
    if (tun == 'user'){
      h1 <- NA
      h2 <- round(object$lambda, 5L)
    }
    if (!tun == 'user'){
      h1 <- round(as.numeric(object$misc$testError), 6L)
      h2 <- round(as.numeric(object$lambda), 5L)
    }
    penalty <- list(penalty = noquote("SERP"),
                    tuneMethod = noquote(tun), value = h1, lambda = h2)
    object$penalty <- penalty
  }
  null <- update(object, ~1)
  object$null.logLik <- null$logLik
  object$null.deviance <- null$deviance
  object$hessian <- H
  object$coefficients <- coefs
  object$expcoefs <- expcoefs
  class(object) <- "summary.serp"
  object
}

#' Prediction from fitted serp model
#'
#' This function takes a fitted \code{serp} object produced by serp() and
#' produces predicted values. Type of predictions returned include response,
#' link and class. Prediction is also possible with new set of values having
#' the same column names as in the original values used for the model fit.
#'
#' @param object An object of class \code{serp}.
#' @param type could be any of these: response, link or terms.
#' @param newdata fresh dataset with all relevant variables.
#' @param ... additional arguments.
#' @return A vector of predicted classes with \code{type} equal to 'class'
#' or a dataframe of predicted values for \code{type} equal to 'response'
#' and 'link'.
#' @seealso \code{\link{anova.serp}}, \code{\link{summary.serp}},
#' \code{\link{confint.serp}}, \code{\link{vcov.serp}}
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "penalize",
#'            reverse = TRUE, link = "logit", tuneMethod = "user",
#'            lambda = 1, data = wine)
#'
#' head(predict(m, type = "link"))
#' head(predict(m, type = "response"))
#' predict(m, type = "class")
#'
#' n.wine <- wine[1:20,]
#' predict(m, newdata = n.wine, type = "class")
#'
#' @export
#'
predict.serp <- function(
  object,
  type = c("link", "response", "class"),
  newdata = NULL,
  ...)
{
  if (!inherits(object, "serp"))
    stop("object must be of class \"serp\"", call. = FALSE)
  type <- match.arg(type)
  nL <- object$ylev
  control <- object$control
  resp <- object$fitted.value
  depvar <- object$model[,1L]
  xpred <- model.matrix(object$Terms, object$model)
  if (!is.factor(depvar)) depvar <- factor(depvar)
  if (!(is.null(newdata))){
    newnames <- colnames(newdata)
    oldnames <- all.vars(object$Terms)
    if (!all(is.element(oldnames, newnames)))
      stop("names in newdata do not match previous names", call. = FALSE)
    m <- newdata[,oldnames]
    xpred <- model.matrix(object$Terms, newdata)

    if (!dim(xpred)[1] == 1L)
      x <- xpred[,-1L]
    else x <- subset(xpred, select = -1)
    vnull <- ifelse((dim(x)[2] == 1L), TRUE, FALSE)
    useout <- TRUE
    nL <- object$ylev
    obs <- nrow(x)
    xlst <- formxL(xpred, nL, object$slope, object$globalEff,
                   object$model, vnull)
    xMat <- do.call(rbind, xlst)
    linkf <- lnkfun(object$link)
    npar <- length(object$coef)
    coln <- colnames(xpred)[-1L]
    cofx <- est.names(object$coef, object$slope, object$globalEff,
                      coln, x, m, npar, xMat, nL, vnull, useout)
    cofx <- if (object$reverse) -object$coef else object$coef
    if (length(object$coef) == ncol(xMat))
      suppressWarnings(
        resp <- prlg(cofx, xMat, obs, yMtx = NULL,
                     penx = NULL, linkf, control = NULL, wt = NULL)$exact.pr
      )
    if (dim(xpred)[1] == 1L){
      resp <- t(data.frame(resp))
      row.names(resp) <- NULL
    }
  }
  switch(
    type,
    response = {pred <- data.frame(resp)
    colnames(pred) <- levels(droplevels(depvar))},
    link = {cms <- t(apply(resp, 1, cumsum))[,-nL]
    if (dim(xpred)[1] == 1L){
      cms <- t(data.frame(cms))
      row.names(cms) <- NULL
    }
    if (!object$reverse){
      pred <- data.frame(log(cms/(1 - cms)))
      colnames(pred) <- sprintf("%slink(P[Y<=%d])",
                                object$link,1:(nL - 1))
    } else {
      pred <- - data.frame(log(cms/(1 - cms)))
      colnames(pred) <- sprintf("%slink(P[Y>=%d])",
                                object$link,2:nL)
    }},
    class = {ylevs <- levels(droplevels(depvar))
    pred <- factor(max.col(resp), levels = seq_along(ylevs),
                   labels = ylevs)})
  pred
}

#' AIC for a fitted serp object
#'
#' Returns the akaike information criterion of a fitted object of class
#' \code{serp}. For the penalized slope, the effective degrees of freedom (edf)
#' is obtained from the trace of the generalized hat matrix which depends on
#' the tuning parameter.
#'
#' @param object An object of class \code{serp}.
#' @param k fixed value equal to 2.
#' @param ... additional arguments.
#' @return A single numeric value of the model AIC.
#' @seealso \code{\link{serp}}, \code{\link{BIC.serp}}, \code{\link{coef.serp}},
#' \code{\link{logLik.serp}},
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "parallel", link = "probit",
#'           data = wine)
#' AIC(m)
#' @export
#'
AIC.serp <- function (object, ..., k=2)
{
  edf <- object$edf
  -2 * object$logLik + k * edf
}

#' BIC for a fitted serp object
#'
#' Returns the bayesian information criterion of a fitted object of class
#' \code{serp}. For the penalized slopes, the effective degrees of freedom (edf)
#' is obtained from the trace of the generalized hat matrix which depends on
#' the tuning parameter.
#'
#' @param object An object of class \code{serp}.
#' @param ... additional arguments.
#' @return A single numeric value of the model.
#' @seealso \code{\link{serp}}, \code{\link{AIC.serp}}, \code{\link{coef.serp}},
#' \code{\link{logLik.serp}},
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "parallel", link = "loglog",
#'           data = wine)
#' BIC(m)
#' @export
#'
BIC.serp <- function (object, ...)
{
  edf <- object$edf
  -2 * object$logLik + log(object$nobs) * edf
}

#' Coefficients for a fitted serp object
#'
#' Returns the coefficients of a fitted object of class \code{serp}.
#'
#' @param object An object of class \code{serp}.
#' @param ... additional arguments.
#' @return A vector of model coefficients.
#' @seealso \code{\link{serp}}, \code{\link{AIC.serp}}, \code{\link{BIC.serp}},
#' \code{\link{logLik.serp}}
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "parallel", link = "loglog",
#'           data = wine)
#' coef(m)
#'
#' @export
#'
coef.serp <- function(object, ...)
{
  object$coef
}
#' @rdname coef.serp
coefficients.serp <- coef.serp

#' Log-likelihood for a fitted serp object
#'
#' Returns the Log-likelihood for a fitted object of class \code{serp}.
#'
#' @param object An object of class \code{serp}.
#' @param ... additional arguments.
#' @return A single numeric value of model log-likelihood
#' @seealso \code{\link{serp}}, \code{\link{AIC.serp}}, \code{\link{BIC.serp}},
#' \code{\link{coef.serp}}
#' @examples
#' library(serp)
#' m <- serp(rating ~ temp + contact, slope = "parallel", link = "loglog",
#'           data = wine)
#' logLik(m)
#' @export
#'
logLik.serp <- function(object, ...)
{
  object$logLik
}
