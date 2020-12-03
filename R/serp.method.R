#' Print method for a serp object.
#'
#' Prints the data frame returned by the \code{summary.serp} method.
#'
#' @param x A serp object.
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
  if (!inherits(x, "serp")) stop("not a \"serp\" object", call. = FALSE)
  cat("call:\n")
  print(x$call)
  cat("\n")
  cat("coefficients:\n")
  print(x$coef)
  cat("\nloglik:", x$logLik, " ","aic:", x$aic, "\n")
  if (x$slope == 'penalize'){
    tun <- x$tuning
    if (tun == 'manual'){
      h1 <- NA
      h2 <- round(x$lambda,5)
    }
    if (!tun == 'manual'){
      h1 <- round(as.numeric(x$value),6)
      h2 <- round(as.numeric(x$lambda),5)
      if (is.na(x$logLik)) h1 <- h2 <- NA
    }
    cat("\n")
    cat("\nPenalization details:")
    cat("\npenalty: SERP","  " ,"tuning:", tun, "  " ,"value:", h1, "  ","lambda:", h2, "\n")
  }
  na.ac <- length(x$na.action)
  if (na.ac > 0){
    cat("\n")
    if (na.ac==1){cat("\n---",na.ac,"observation deleted due to missingness","---", "\n")
    } else {
      cat("\n")
      cat("\n---",na.ac,"observations deleted due to missingness","---", "\n")
    }
  }
}

#' Summary method for a serp object.
#'
#' Summarizes the results of the fitted model in a dataframe.
#'
#' @param object A "serp" object.
#' @param ... Not used. Additional summary arguments.
#' @return A data frame of the fitted model summary.
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#' @export
#'
summary.serp <- function(object, ...){
  if (!inherits(object, "serp")) stop("not a \"serp\" object", call. = FALSE)
  coef <- object$coef
  nL <- object$ylev
  coefs <- matrix(0, length(coef), 4L, dimnames = list(names(coef), NULL))
  coefs[, 1L] <- coef
  H <- cbind(object$hess[,c(1:ncol(object$hess))])
  dimnames(H) <- list(names(object$coef), names(object$coef))
  cholHx <- try(chol(H), silent = TRUE)
  if (inherits(cholHx, "try-error"))
    cholHx <- covx <- diag(NA, dim(H))
  else
    covx   <- chol2inv(cholHx)
  coefs[, 2L] <- se.est <- sqrt(diag(covx))
  coefs[, 3L] <- z.value <- coef/se.est
  coefs[, 4L] <- pvalue  <- 2 * pnorm(abs(z.value), lower.tail = FALSE)
  colnames(coefs) <- c("Estimate", "Std Error", "z value", "Pr(>|z|)")
  expcoefs <- exp(coef(object)[-c(1:(nL-1))])
  if (object$slope == 'penalize'){
    tun <- object$tuning
    if (tun == 'manual'){
      h1 <- NA
      h2 <- round(object$lambda, 5L)
    }
    if (!tun == 'manual'){
      h1 <- round(as.numeric(object$value), 6L)
      h2 <- round(as.numeric(object$lambda), 5L)
      if (is.na(object$logLik)) h1 <- h2 <- NA
    }
    penalty <- list(penalty = noquote("SERP"), tuning = noquote(tun), value = h1, lambda = h2)
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
  if (!inherits(x, "summary.serp")) stop("not a \"serp\" object", call. = FALSE)
  object <- x
  cat("call:\n")
  print(object$call)
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
  if (object$slope == 'penalize'){
    tun <- object$tuning
    if (tun =='manual'){
      h1 <- NA
      h2 <- round(object$lambda, 5L)
    }
    if (!tun == 'manual'){
      h1 <- round(as.numeric(object$value), 6L)
      h2 <- round(as.numeric(object$lambda), 5L)
      if (is.na(object$logLik)) h1 <- h2 <- NA
    }
    cat("\n")
    cat("\nPenalization details:")
    cat("\npenalty: SERP","  " ,"tuning:", tun, "  " ,"value:", h1, "  ","lambda:", h2, "\n")
  }
  if (na.ac > 0){
    cat("\n")
    if (na.ac == 1){
      cat("\n---",na.ac,"observation deleted due to missingness","---")
    } else{
      cat("\n")
      cat("\n---",na.ac,"observations deleted due to missingness","---")
    }
  }
  invisible(object)
}

#' Predict method for a serp object.
#'
#' Returns the predicted probabilities, link and class for an object of class "serp".
#'
#' @param object A "serp" object.
#' @param type any of response, link or terms.
#' @param newdata fresh dataset with all relivant variables.
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
  if (!inherits(object, "serp")) stop("not a \"serp\" object", call. = FALSE)
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
      stop("variables in new data different from those in the model")
    m <- newdata[,oldnames]
    xpred <- model.matrix(object$Terms, newdata)
    if (!dim(xpred)[1] == 1L) x <- xpred[,-1L] else x <- subset(xpred, select = -1)
    vnull <- ifelse((dim(x)[2] == 1L), TRUE, FALSE)
    useout <- TRUE
    ym <- oldnames[1L]
    y <- newdata[,ym]
    if (!is.factor(y)) y <- factor(y)
    wt <- rep(1, length(y))
    if (!is.factor(y))
      stop("response must be a factor", call. = FALSE)
    nL <- length(levels(y))
    if (nL <= 2)
      stop("response must have more than two levels", call. = FALSE)
    obs <- length(y)
    if (obs != dim(x)[1])
      stop("variable lengths unequal", call. = FALSE)
    yMtx <- yMx(y, obs, nL)
    if (object$ylev != nL){
      yint <- intersect(levels(y),levels(droplevels(depvar)))
      yMtx <- yMtx[ ,yint]
      nL <- ncol(yMtx)
    }
    if (nL > 2L){
      xlst <- formxL(xpred, nL, object$slope, object$globalvar, object$model, vnull)
      xMat <- do.call(rbind, xlst)
      linkf <- lnkfun(object$link)
      npar <- length(object$coef)
      coln <- colnames(xpred)[-1L]
      cofx <- est.names(object$coef, object$slope, object$globalvar,
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
    type, response= {pred <- data.frame(resp)
    colnames(pred) <- levels(droplevels(depvar))},
    link = {cms <- t(apply(resp, 1, cumsum))[,-nL]
    if (dim(xpred)[1] == 1L){
      cms <- t(data.frame(cms))
      row.names(cms) <- NULL
    }
    pred <- data.frame(log(cms/(1-cms)))
    colnames(pred) <- sprintf("%slink(P[Y<=%d])", object$link,1:(nL-1))},
    class = {ylevs <- levels(droplevels(depvar))
    pred <- factor(max.col(resp), levels = seq_along(ylevs), labels = ylevs)}
  )
  pred
}

#' AIC for a serp object.
#'
#' @param object A "serp" object.
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
#' @param object object of class 'serp'.
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
#' @param object A "serp" object.
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
#' @param object A "serp" object.
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

#' Performance metrics for categorical models
#'
#' @description Calculates performance metrics of fitted categorical models,
#' including binary and multi-categorical models.
#'
#' @usage errorMetrics(actual, predicted, model= c("multiclass", "binary"),
#'                     type= c("brier", "logloss", "misclass"), eps=1e-15)
#' @param actual vector of actual values observed
#' @param predicted predicted probability matrix of a categorical model or a vector of fitted values for binary models.
#' @param model specifies whether multi-categorical or binary model
#' @param type specifies type of error metrics
#' @param eps a near-zero value introduced only if the fitted probabilities goes beyond
#' specified threshold. Helps to minimize the chances of running into
#' numerical problems
#' @return A numeric value of computed performance metric determining how
#' good a categorical model is compare to competing models.
#' \describe{
#'   \item{brier}{the brier score of fitted model.}
#'   \item{logloss}{the logloss of fitted model.}
#'   \item{misclass}{the misclassification error of fitted model.}
#'}
#'
#' @seealso
#' \code{\link{serp}}
#'
#' @examples
#' \dontrun{
#' f1 <- serp(rating ~ temp + contact, globalvar = ~ temp,
#' slope = "partial", reverse = T, link = "logit",
#' data = wine)
#' y1 <- wine$rating
#' p1 <- predict(f1, type="response")
#'
#' errorMetrics(f1, type = "brier")
#' errorMetrics(f1, type = "logloss")
#' errorMetrics(f1, type = "misclass")
#'
#' ## For other class of model 'actual'  and predicted must be provided
#' set.seed(1)
#' y <- as.factor(rbinom(50,1,0.5))
#' xx <- runif(50)
#' f2 <- glm(y ~ xx, family = binomial(link="logit"))
#' p2 <- f2$fitted.values
#'
#' errorMetrics(actual=y, predicted=p2, model= "binary", type = "brier")
#' errorMetrics(actual=y, predicted=p2, model= "binary", type = "logloss")
#' errorMetrics(actual=y, predicted=p2, model= "binary", type = "misclass")
#' }
#'
#' @export
#'
errorMetrics <- function(actual, predicted, model= c("multiclass", "binary"),
                         type= c("brier", "logloss", "misclass"), eps=1e-15)
{
  model <- match.arg(model)
  type <- match.arg(type)
  y <- actual

  if (inherits(y, "serp")){
    obj <- y
    y <- factor(y$model[,1L])
    pred_y <- obj$fitted.values
  }else{
    pred_y <- predicted
  }
  if(!is.factor(y)) stop("y must be a factor")
  y <- droplevels(y)
  nL <- nlevels(y)
  obs <- length(y)
  NAs <- 0
  if(nL < 2)stop("y must have two or more levels")
  if(model=="multiclass"){
    if(is.null(ncol(pred_y)))
      stop("supply either a matrix or dataframe of fitted values")
    if(nrow(pred_y) != obs || max(unclass(y)) != ncol(pred_y) || min(unclass(y)) < 0)
      stop("y levels not equal to the number of columns of fitted values, or unequal lengths of observations")
  }else{
    if(!is.null(dim(pred_y)))
      stop("supply a vector of fitted values")
    if(!is.numeric(pred_y)) stop("supply a numeric vector of fitted values")
    if(length(pred_y) != obs) stop("lengths of y and fitted values unequal")
  }
  yprob <- cbind(y, pred_y)
  if(any(is.na(pred_y))){
    yprob <- na.omit(yprob)
    y <- yprob[,1L]
    pred_y <- yprob[,-1L]
    newobs <- nrow(yprob)
    NAs <- obs - newobs
    obs <- newobs
  }
  pred_y[pred_y > 1- eps] <- 1 - eps
  pred_y[pred_y < eps] <- eps
  switch(
    model,
    binary={
      if(type=="brier"){
        y <- as.integer(y) - 1L
        error <- sum((y - pred_y)^{2})/obs
      }
      if(type=="logloss"){
        y <- as.integer(y) - 1L
        error <- - (sum(y * log(pred_y) + (1 - y) * log(1 - pred_y))) / length(y)
      }
      if(type=="misclass"){
        y <- factor(y)
        ylevs <- levels(y)
        rr <- factor(max.col(pred_y), levels = seq_along(ylevs), labels = ylevs)
        error <- mean(y!=rr)
      }
    },
    multiclass={
      ym <- matrix(0, nrow=obs, ncol=nL, dimnames=list(NULL, levels(y)))
      yi <- as.integer(y)
      ym[cbind(1:obs, yi)] <- 1

      if(type=="brier"){
        rs <- rowSums(ym)
        error <- sum(ym * (1 - pred_y)^2 + (rs - ym) * pred_y^2) / sum(rs)
      }
      if(type=="logloss"){
        error <- -sum(ym * log(pred_y))/nrow(pred_y)
      }
      if(type=="misclass"){
        rr <- apply(pred_y, 1, which.max)
        hh <- sapply(1:nrow(ym), function(i) sum(ym[i, -rr[i]]))
        error <- sum(hh) / sum(ym)
      }
    })
  if(NAs > 0) cat("\n---",NAs,"observation(s) deleted due to missingness","---", "\n")
  error
}

#' Control parameters for serp fit.
#'
#' Default control parameters for serp fit. User supplied control parameters could be
#' specified in the main function.
#'
#' @usage serp.control(maxits = 500, eps = 1e-06, maxpen = 1e5, nFold = 5)
#' @param maxits the maximum number of iteration.
#' @param eps threshold value during optimization at which the iteration routine terminates. In other words, when the reported change in the log-likelihood goes below this threshold, convergence is achieved.
#' @param maxpen the upper end point of the interval from zero to be searched for a minimum tuning parameter.
#' @param nFold the number of k-fold cross validation for the "cv" \code{tuning} method, with default value k = 5.
#' @return NULL
#' @seealso
#' \code{\link{serp}}
#' @examples
#' # See serp() documentation for examples.
#'
#' @export
#'
serp.control <- function (maxits = 5e02, eps = 1e-06, maxpen = 1e05,
                          nFold = 5e0)
{
  if (!is.numeric(eps) || eps <= 0)
    stop("value of eps must be > 0")
  if (!is.numeric(maxits) || maxits <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(maxpen) || maxpen <= 0)
    stop("value of maxpen must be > 0")
  if (is.null(maxits)) maxits <- 5e02
  if (is.null(eps)) eps <- 1e-06
  if (is.null(maxpen)) maxpen <- 1e05
  if(!(is.numeric(nFold)) | (nFold < 2L) | (nFold > 10L)|(round(nFold) != nFold))
    stop("nFold should be numeric, non-decimal, and between 2 and 10 inclusive.", call. = FALSE)
  if (is.null(nFold)) nFold <- 5e0
  list(maxits = maxits, eps = eps, maxpen = maxpen, nFold = nFold)
}
