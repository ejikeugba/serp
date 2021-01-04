#' Print method for a serp object.
#'
#' Prints the data frame returned by the \code{summary.serp} method.
#'
#' @param x An object of class serp.
#' @param ... additional arguments.
#' @return NULL
#' @seealso
#' \code{\link{serp}}
#'
#' @examples
#' # See serp() documentation for examples.
#'
#' @export
#'
print.serp <- function(x, ...)
{
  if (!inherits(x, "serp"))
    stop("not a \"serp\" object", call. = FALSE)
  object <- x
  cat("\ncall:\n")
  print(object$call)
  if (object$slope == 'penalize'){
    tun <- object$tuning
    if (tun =='manual'){
      h1 <- " "
      h2 <- round(object$lambda, 5L)
    }
    if (!tun == 'manual'){
      h1 <- round(as.numeric(object$value), 6L)
      h2 <- round(as.numeric(object$lambda), 5L)
      if (is.na(object$logLik)) h1 <- h2 <- NA
    }
    if (object$tuning == "cv" && object$cverror == "brier")
      cverror <- "cv (brier)"
    if (object$tuning == "cv" && object$cverror == "logloss")
      cverror <- "cv (logloss)"
    cat("\nPenalization:")
    if (!length(colnames(object$model)[-1L]) ==
        length(object$global)){
      cat("\npenalty:", "  "," SERP")
      cat("\ntuning:","    " , tun)
      if (object$tuning == "cv")
        cat("\ncverror:",  "   ", object$cverror)
      cat("\nvalue:","     " , h1)
      cat("\nlambda:","    " , h2)
    } else cat("\nwith no subject-specific effects found, ",
               "parallel slope(s) returned\n")
    cat("\n")
  }
  cat("\ncoefficients:\n")
  print(object$coef)
  cat("\nloglik:", object$logLik, " ","aic:", object$aic, "\n")
  na.ac <- length(object$na.action)
  if (na.ac > 0){
    cat("\n")
    if (na.ac==1)
    {
      cat("\n---",na.ac,"observation deleted due to ",
          "missingness","---", "\n")
    } else {
      cat("\n")
      cat("\n---",na.ac,"observations deleted due to ",
          "missingness","---", "\n")
    }
  }
}

#' Summary method for a serp object.
#'
#' Summarizes the results of the fitted model in a dataframe.
#'
#' @param object An object of class serp.
#' @param ... Not used. Additional summary arguments.
#' @return A data frame of the fitted model summary.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#' @export
#'
summary.serp <- function(object, ...){
  if (!inherits(object, "serp"))
    stop("not a \"serp\" object", call. = FALSE)
  coef <- object$coef
  nL <- object$ylev
  coefs <- matrix(0, length(coef), 4L,
                  dimnames = list(names(coef), NULL))
  coefs[, 1L] <- coef
  H <- cbind(object$hess[,c(1:ncol(object$hess))])
  dimnames(H) <- list(names(object$coef), names(object$coef))
  cholHx <- try(chol(H), silent = TRUE)
  if (inherits(cholHx, "try-error"))
    cholHx <- covx <- diag(NA, dim(H))
  else covx   <- chol2inv(cholHx)
  coefs[, 2L] <- se.est <- sqrt(diag(covx))
  coefs[, 3L] <- z.value <- coef/se.est
  coefs[, 4L] <- pvalue  <- 2 * pnorm(abs(z.value), lower.tail = FALSE)
  colnames(coefs) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)")
  expcoefs <- exp(coef(object)[-c(1:(nL-1))])
  if (object$slope == 'penalize'){
    tun <- object$tuning
    if (object$tuning == "cv" && object$cverror == "brier")
      tun <- "cv (brier)"
    if (object$tuning == "cv" && object$cverror == "logloss")
      tun <- "cv (logloss)"
    if (tun == 'manual'){
      h1 <- NA
      h2 <- round(object$lambda, 5L)
    }
    if (!tun == 'manual'){
      h1 <- round(as.numeric(object$value), 6L)
      h2 <- round(as.numeric(object$lambda), 5L)
      if (is.na(object$logLik)) h1 <- h2 <- NA
    }
    penalty <- list(penalty = noquote("SERP"),
                    tuning = noquote(tun), value = h1, lambda = h2)
    object$penalty <- penalty
  }
  object$hess <- H
  object$coef <- coefs
  object$expcoefs <- expcoefs
  class(object) <- "summary.serp"
  object
}

#' @export
#'
print.summary.serp <- function(x, ...){
  if (!inherits(x, "summary.serp"))
    stop("not a \"serp\" object", call. = FALSE)
  object <- x
  cat("\ncall:\n")
  print(object$call)
  if (object$slope == 'penalize'){
    tun <- object$tuning
    if (tun =='manual'){
      h1 <- " "
      h2 <- round(object$lambda, 5L)
    }
    if (!tun == 'manual'){
      h1 <- round(as.numeric(object$value), 6L)
      h2 <- round(as.numeric(object$lambda), 5L)
      if (is.na(object$logLik)) h1 <- h2 <- NA
    }
    if (object$tuning == "cv" && object$cverror == "brier")
      cverror <- "cv (brier)"
    if (object$tuning == "cv" && object$cverror == "logloss")
      cverror <- "cv (logloss)"
    cat("\nPenalization:")
    if (!length(colnames(object$model)[-1L])==length(object$global)){
      cat("\npenalty:", "  "," SERP")
      cat("\ntuning:","    " , tun)
      if (object$tuning == "cv") cat("\ncverror:",  "   ", object$cverror)
      cat("\nvalue:","     " , h1)
      cat("\nlambda:","    " , h2)
    } else cat("\nwith no subject-specific effects found, ",
               "parallel slope(s) returned\n")
    cat("\n")
  }
  nL <- object$ylev
  coef <- as.data.frame(object$coef)
  df <- (object$nobs*(object$ylev-1)) - length(coef[,1L])
  na.ac <- length(object$na.action)
  cat("\nCoefficients:\n")
  printCoefmat(coef, digits = 4L, signif.stars = TRUE,
               na.print = "NA", P.values = TRUE, has.Pvalue = TRUE,
               signif.legend = TRUE)
  cat("\nNumber of iterations:", object$iter, "\n")
  cat("\nLoglik:", object$logLik, "on", df, "degrees of freedom","\n")
  cat("\nAIC:", object$aic)
  cat("\n")
  if (!ncol(object$model) == 1L){
    cat("\nExponentiated coefficients:\n")
    print(object$expcoefs)
  }
  if (na.ac > 0){
    cat("\n")
    if (na.ac == 1){
      cat("\n---",na.ac,"observation deleted ",
          "due to missingness","---")
    } else{
      cat("\n")
      cat("\n---",na.ac,"observations deleted ",
          "due to missingness","---")
    }
  }
  invisible(object)
}

#' Predict method for a serp object.
#'
#' Returns the predicted probabilities, link and class
#' for an object of class "serp".
#'
#' @param object An object of class serp.
#' @param type any of response, link or terms.
#' @param newdata fresh dataset with all relevant variables.
#' @param ... additional arguments.
#' @return The object returned depends on \code{type}.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#' @export
#'
predict.serp <- function(object,
                         type = c("link", "response", "class"),
                         newdata=NULL, ...)
{
  if (!inherits(object, "serp"))
    stop("not a \"serp\" object", call. = FALSE)
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
    if (all(!is.element(oldnames, newnames)))
      stop("variables in new data different ",
           "from those in the model")
    m <- newdata[,oldnames]
    xpred <- model.matrix(object$Terms, newdata)
    if (!dim(xpred)[1] == 1L)
      x <- xpred[,-1L]
    else x <- subset(xpred, select = -1)
    vnull <- ifelse((dim(x)[2] == 1L), TRUE, FALSE)
    useout <- TRUE
    ym <- oldnames[1L]
    y <- newdata[,ym]
    if (!is.factor(y)) y <- factor(y)
    wt <- rep(1, length(y))
    if (!is.factor(y))
      stop("response must be a factor",
           call. = FALSE)
    nL <- length(levels(y))
    if (nL <= 2)
      stop("response must have more than two levels",
           call. = FALSE)
    obs <- length(y)
    if (obs != dim(x)[1])
      stop("variable lengths unequal",
           call. = FALSE)
    yMtx <- yMx(y, obs, nL)
    if (object$ylev != nL){
      yint <- intersect(levels(y),levels(droplevels(depvar)))
      yMtx <- yMtx[ ,yint]
      nL <- ncol(yMtx)
    }
    if (nL > 2L){
      xlst <- formxL(xpred, nL, object$slope, object$global,
                     object$model, vnull)
      xMat <- do.call(rbind, xlst)
      linkf <- lnkfun(object$link)
      npar <- length(object$coef)
      coln <- colnames(xpred)[-1L]
      cofx <- est.names(object$coef, object$slope, object$global,
                        coln, x, m, npar, xMat, nL, vnull, useout)
      if (length(object$coef) == ncol(xMat))
        suppressWarnings(
          resp <- prlg(cofx, xMat, obs, yMtx=NULL,
                       penx=NULL, linkf, control=NULL, wt=NULL)$prob
        )
      if (dim(xpred)[1] == 1L){
        resp <- t(data.frame(resp))
        row.names(resp) <- NULL
      }
    }
  }
  switch(
    type,
    response= {pred <- data.frame(resp)
    colnames(pred) <- levels(droplevels(depvar))},
    link = {cms <- t(apply(resp, 1, cumsum))[,-nL]
    if (dim(xpred)[1] == 1L){
      cms <- t(data.frame(cms))
      row.names(cms) <- NULL}
    pred <- data.frame(log(cms/(1-cms)))
    colnames(pred) <- sprintf("%slink(P[Y<=%d])",
                              object$link,1:(nL-1))},
    class = {ylevs <- levels(droplevels(depvar))
    pred <- factor(max.col(resp), levels = seq_along(ylevs),
                   labels = ylevs)})
  pred
}

#' AIC for a serp object.
#'
#' @param object An object of class serp.
#' @param k fixed value equal to 2.
#' @param ... additional arguments.
#' @return A single value of model AIC.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#' @export
#'
AIC.serp <- function (object, ..., k=2)
{
  edf <- object$edf
  -2 * object$logLik + k * edf
}

#' BIC for a serp object.
#'
#' @param object An object of class serp.
#' @param ... additional arguments.
#' @return A single value of model BIC.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#' @export
#'
BIC.serp <- function (object, ...)
{
  edf <- object$edf
  -2 * object$logLik + log(object$nobs) * edf
}

#' Coefficients for a serp object.
#'
#' @param object An object of class serp.
#' @param ... additional arguments.
#' @return A vector of model coefficients.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#' @export
#'
coef.serp <- function(object, ...)
{
  object$coef
}

#' Log-likelihood for a serp object.
#'
#' @param object An object of class serp.
#' @param ... additional arguments.
#' @return A single value of model log-likelihood
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#' @export
#'
logLik.serp <- function(object, ...)
{
  object$logLik
}
