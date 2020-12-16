#' Performance metrics for categorical models
#'
#' @description Calculates performance metrics of fitted categorical models,
#' including binary and multi-categorical models.
#'
#' @usage errorMetrics(actual, predicted, model= c("multiclass", "binary"),
#'                     type= c("brier", "logloss", "misclass"),
#'                     eps=.Machine$double.eps)
#' @param actual vector of actual values observed
#' @param predicted predicted probability matrix of a categorical model or a vector of fitted values for binary models.
#' @param model specifies whether multi-categorical or binary model
#' @param type specifies type of error metrics
#' @param eps a near-zero value introduced only if the fitted probabilities goes beyond
#' specified threshold. Helps to minimize the chances of running into
#' numerical problems.
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
#' ## For non-serp object of class, 'actual'  and 'predicted' values
#' ## must be provided
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
errorMetrics <- function(actual, predicted,
                         model = c("multiclass", "binary"),
                         type = c("brier", "logloss", "misclass"),
                         eps = .Machine$double.eps)
{
  model <- match.arg(model)
  type <- match.arg(type)
  y <- actual
  if (inherits(y, "serp")){
    obj <- y
    y <- factor(y$model[,1L])
    pred_y <- obj$fitted.values
  }
  else
    pred_y <- predicted
  if (!is.factor(y))
    stop("'actual' must be a factor")
  y <- droplevels(y)
  nL <- nlevels(y)
  obs <- length(y)
  NAs <- 0
  if (nL < 2)
    stop("actual observations must have two or more levels")
  if (model=="multiclass"){
    if (is.null(ncol(pred_y)))
      stop("supply either a matrix or ",
           "dataframe of fitted values")
    if (nrow(pred_y) != obs || max(unclass(y)) != ncol(pred_y)
        || min(unclass(y)) < 0)
      stop("levels actual observations not equal to the number of ",
           "columns of fitted values, or unequal ",
           "lengths of observations")
  }else{
    if (!is.null(dim(pred_y)))
      stop("supply a vector of fitted values")
    if (!is.numeric(pred_y))
      stop("supply a numeric vector of fitted values")
    if (length(pred_y) != obs)
      stop("lengths of actual observations and fitted values unequal")
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
  pred_y[pred_y < eps] <- eps
  switch(
    model,
    binary={
      if (type=="brier"){
        y <- as.integer(y) - 1L
        error <- sum((y - pred_y)^{2})/obs}
      if (type=="logloss"){
        y <- as.integer(y) - 1L
        error <- - (sum(y * log(pred_y) + (1 - y) *
                          log(1 - pred_y))) / length(y)}
      if (type=="misclass"){
        y <- factor(y)
        ylevs <- levels(y)
        rr <- factor(max.col(pred_y), levels = seq_along(ylevs),
                     labels = ylevs)
        error <- mean(y!=rr)}
    },
    multiclass={
      ym <- matrix(0, nrow=obs, ncol=nL,
                   dimnames=list(NULL, levels(y)))
      yi <- as.integer(y)
      ym[cbind(1:obs, yi)] <- 1
      if (type=="brier"){
        rs <- rowSums(ym)
        error <- sum(ym * (1 - pred_y)^2 + (rs - ym) *
                       pred_y^2) / sum(rs)}
      if (type=="logloss"){
        error <- -sum(ym * log(pred_y))/nrow(pred_y)}
      if (type=="misclass"){
        rr <- apply(pred_y, 1, which.max)
        hh <- sapply(1:nrow(ym), function(i) sum(ym[i, -rr[i]]))
        error <- sum(hh) / sum(ym)}
    })
  if (NAs > 0)
    cat("\n---",NAs,"observation(s) deleted ",
        "due to missingness","---", "\n")
  error
}
