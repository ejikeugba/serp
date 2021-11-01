#' Performance metrics for categorical response models
#'
#' @description Calculates the performance metrics of fitted binary and
#' multi-categorical response models. Available measures include: brier score,
#' logloss and misclassification error.
#'
#' @usage errorMetrics(
#'              actual,
#'              predicted,
#'              type= c("brier", "logloss", "misclass"),
#'              control = list())
#' @param actual vector of actual values observed
#' @param predicted predicted probability matrix of a categorical model or a
#'   vector of fitted values for binary models.
#' @param type specifies the type of error metrics.
#' @param control A list of control parameters to replace default values
#'   returned by \code{serp.control}. 'misclass.thresh' resets the default
#'   misclassification error threshold, while 'minP' assigns a near-zero constant
#'   value to the predicted values beyond certain threshold, to forestall chances
#'   of numerical problems.
#' @details A numeric value of the selected error type determining how
#' good a categorical model compares to competing models.
#' @return \item{value}{the value of error measure computed.}
#' @return \item{type}{the error measure used: any of brier, logLoss or misclassification error.}
#' @return \item{threshold}{the misclassification threshold.}
#'
#' @seealso
#' \code{\link{serp}}, \code{\link{anova.serp}}, \code{\link{confint.serp}},
#' \code{\link{vcov.serp}}
#'
#' @examples
#'
#' require(serp)
#'
#' m1 <- serp(rating ~ temp + contact, slope = "parallel", link = "logit", data = wine)
#' errorMetrics(m1, type = "brier")
#'
#' ## objects of class other than \code{serp} require the actual
#' ## observations with corresponding predicted values supplied.
#'
#' set.seed(2)
#' n <- 50
#' y <- as.factor(rbinom(n, 1, 0.5))
#' m2 <- glm(y ~ rnorm(n), family = binomial())
#' ft <- m2$fitted.values
#'
#' errorMetrics(actual=y, predicted=ft, type = "logloss")
#' errorMetrics(actual=y, predicted=ft, type = "misclass")
#'
#' ## Reset classification threshold
#' errorMetrics(actual=y, predicted=ft, type = "misclass",
#' control = list(misclass.thresh=0.4))
#'
#' @export
#'
errorMetrics <- function(
  actual,
  predicted,
  type = c("brier", "logloss", "misclass"),
  control = list())
{
  type <- match.arg(type)
  control <- do.call("serp.control", control)
  y <- actual
  if (inherits(y, "serp")){
    obj <- y
    y <- factor(y$model[,1L])
    pred_y <- obj$fitted.values
  } else {
    if (missing(actual) || missing(predicted))
      stop("please provide actual and predicted values of fit!")
    y <- actual
    pred_y <- predicted
  }
  if (!is.factor(y))
    stop("'actual' must be a factor")
  y <- droplevels(y)
  nL <- nlevels(y)
  obs <- length(y)
  NAs <- 0
  if (nL < 2)
    stop("actual observations must have two or more levels")
  category <- if (nL > 2) "multiclass" else "binary"
  if (category=="multiclass"){
    if (is.null(ncol(pred_y)))
      stop("supply either a matrix or ",
           "dataframe of fitted values")
    if (nrow(pred_y) != obs || max(unclass(y)) != ncol(pred_y)
        || min(unclass(y)) < 0)
      stop("levels of actual observations not equal to the number of ",
           "columns of fitted values, or unequal ",
           "lengths of observations")
  }else{
    if (!is.null(dim(pred_y)))
      stop("supply a vector of fitted values")
    if (!is.numeric(pred_y))
      stop("supply a numeric vector of fitted values")
    if (length(pred_y) != obs)
      stop("unequal lengths of actual and fitted values")
  }
  yprob <- cbind(y, pred_y)
  if (any(is.na(pred_y))){
    yprob <- na.omit(yprob)
    y <- yprob[,1L]
    pred_y <- yprob[,-1L]
    newobs <- nrow(yprob)
    NAs <- obs - newobs
    obs <- newobs
  }
  eps <- control$minP
  pred_y[pred_y < eps] <- eps
  pred_y[pred_y > 1-eps] <- 1-eps
  pred_y <- if (category=="multiclass") pred_y/rowSums(pred_y) else pred_y
  switch(
    category,
    binary={
      if (type=="brier"){
        y <- as.integer(y) - 1L
        value <- sum((y - pred_y)^{2})/obs}
      if (type=="logloss"){
        y <- as.integer(y) - 1L
        value <- - (sum(y * log(pred_y) + (1 - y) *
                          log(1 - pred_y))) / length(y)}
      if (type=="misclass"){
        y <- factor(y)
        ylevs <- levels(y)
        rr <- ifelse(pred_y > control$misclass.thresh, 1, 0)
        value <- mean(y!=rr)}
    },
    multiclass={
      ym <- matrix(0, nrow=obs, ncol=nL,
                   dimnames=list(NULL, levels(y)))
      yi <- as.integer(y)
      ym[cbind(1:obs, yi)] <- 1
      if (type=="brier"){
        rs <- rowSums(ym)
        value <- sum(ym * (1 - pred_y)^2 + (rs - ym) *
                       pred_y^2) / sum(rs)}
      if (type=="logloss"){
        value <- -sum(ym * log(pred_y))/nrow(pred_y)}
      if (type=="misclass"){
        rr <- apply(pred_y, 1, which.max)
        rr <- rr - min(rr) + 1L
        hh <- vapply(seq_len(nrow(ym)), function(i) sum(ym[i, -rr[i]]),
                     numeric(1))
        value <- sum(hh) / sum(ym)}
    })
  if (NAs > 0)
    warning(NAs," ", "observation(s) deleted ","due to missingness")

  res <- list(value=value, type=type)
  if (category == "binary" && type == "misclass")
    res$threshold <- control$misclass.thresh
  class(res) <- "errorMetrics"
  res
}
