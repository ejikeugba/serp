#' Smooth Effects on Response Penalty for CLM
#'
#' @description Fits cumulative link model (CLM) with smooth-effect-on-response penalty (SERP)
#' via Newton Raphson algorithm. SERP enables the regularization of the parameter space
#' between the general cumulative model and the restricted cumulative model,
#' resulting in all subject-specific effects shrinked to global effects.
#'
#' @usage serp(formula, link = c("logit", "probit","loglog", "cloglog", "cauchit"),
#'             slope = c("penalize", "parallel", "unparallel", "partial"),
#'             tuning = c("deviance","cv","manual","finite"), reverse = FALSE,
#'             lambdagrid=NULL, globalvar = NULL, data, subset, weights=NULL,
#'             weight.type = c("analytic", "frequency"), na.action=NULL,
#'             lambda = NULL, contrasts = NULL, control = list(), ...)
#' @param formula regression formula of the form: response ~ predictors. The response should be a factor (ordered).
#' @param link specifies the link function for the cumulative link model including: logit, probit, complementary log-log, cloglog, cauchit.
#' @param slope "penalize" for penalized unparallel coefficient terms, "unparallel" , "parallel" and "partial" for unpenalized non-parallel, "parallel" and "partial" coefficient terms respectively.
#' @param tuning specifies the method of choosing an optimal shrinkage parameter, including: deviance, cv, manual and finite. i.e., the lambda value along parameter shrinkage path at which the residual deviance of model is minimal or at which the cross-validated prediction error (brier score) is minimal, or alternatively, a user supplied lambda. The "finite" tuning is used to locate the lambda value for which the fit's log-Likelihood is finite.
#' @param reverse false by default, when true the sign of the linear predictor is reversed.
#' @param lambdagrid optional user supplied lambda grid for cv and deviance tuning methods. Negative range of values are not allowed, instead (0, Inf). With large number of predictors and cases in the model, iterations run faster over a short grid interval, i.e., the shorter the grid length the faster iteration.
#' @param globalvar specifies variables to be assigned global effects when \code{slope} is set to \code{partial}. Variables are specified as a formula with an empty left hand side, for instance, globalvar = ~predictors.
#' @param data optional data frame explaining the variables used in the formula.
#' @param subset specifies which subset of the rows of the data should be used for in fit. All observations are used by default.
#' @param weights optional case weights in fitting. Negative weights are not allowed. Defaults to 1.
#' @param weight.type distinguishes between analytic and frequency weights with the former set as default. The latter should be used when weights are mere case counts used to compress the data set.
#' @param na.action a function to filter missing data.
#' @param lambda a user-specified single numeric value for the tuning parameter when the \code{tuning} method is set to manual. A negative value is not allowed.
#' @param contrasts a list of contrasts to be used for some or all of the factors appearing as variables in the model formula.
#' @param control A list of fit control parameters to replace default values returned by \code{serp.control}. Values not set assume default values.
#' @param ... additional arguments.
#' @import stats
#' @importFrom ordinal pgumbel
#' @importFrom ordinal dgumbel
#' @importFrom ordinal qgumbel
#' @seealso \code{\link{anova.serp}}, \code{\link{summary.serp}},
#' \code{\link{predict.serp}}, \code{\link{confint.serp}},
#' \code{\link{vcov.serp}}, \code{\link{errorMetrics}}
#' @details
#' The \code{serp} function fits the cumulative link model (CLM)
#' with smooth-effect-on-response penalty (SERP). The cumulative
#' model developed by McCullagh (1980) is probably most frequently
#' used ordinal model. When motivated by an underlying latent
#' variable, a simple form of the model is expressed as follows:
#'
#' \deqn{P(Y\leq r|x) = F(\delta_{0r} + x^T\delta)}
#'
#' where \eqn{x} is a vector of covariates, \eqn{\delta} a vector
#' of regression parameters and  \eqn{F} a continuous distribution
#' function. The choice of \eqn{F}, of course, results in different
#' forms of the cumulative model. As observed, this assumes that the
#' effect of \eqn{x} does not depend on the  category. However,
#' relaxing such an assumption allows the effect of the covarites to
#' vary across categories, yielding the following general cumulative
#' model:
#'
#' \deqn{P(Y\leq r|x) = F(\delta_{0r} + x^T\delta_{r}),   r=1,\dots,k-1}
#'
#' This model has the property of stochastic ordering (McCullagh, 1980),
#' which implies that  \eqn{\delta_{0r} + x^T\delta_{r-1} < \delta_{0r} + x^T\delta_{r}}
#' holds for all \eqn{x} and all categories  \eqn{r}, since
#' \eqn{P(Y\leq r-1|x) < P(Y\leq r|x)} must hold for all categories.
#' A mixture of the two effect types can also be made, for instance,
#' with the cumulative logit model, such results in the so-called
#' partial proportional odds model. Assuming  \eqn{x} is partitioned
#' into \eqn{x_{1}} and  \eqn{x_{2}}, for all possible link functions
#' this may be expressed as follows:
#'
#' \deqn{P(Y\leq r|x) = F(\delta_{0r} + x^T_{1}\delta + x^T_{2}\delta_{r}), r=1,\dots,k-1}
#'
#' However, with the stochastic ordering constraint on the general
#' cumulative model very often problematic, resulting in unstable
#' likelihoods with ill-conditioned parameter space during the
#' iterative procedure, SERP offers a means of arriving at stable
#' estimates of the general model. It provides a form of regularization
#' that is based on minimizing the penalized log-likelihood:
#'
#' \deqn{l_{p}(\delta)=l(\delta)-J_{\lambda}(\delta)}
#'
#' where  \eqn{l(\delta)}, is the log-likelihood of the general
#' cumulative model and  \eqn{J_{\lambda}(\delta)=\lambda J(\delta)}
#' the penalty function weighted by the turning parameter \eqn{\lambda}.
#' Assuming an ordered categorical outcome  \eqn{Y \in \{1,\dots,k\}},
#' and considering that the corresponding parameters  \eqn{\delta_{1j},\dots \delta_{k-1,j}}
#' vary smoothly over the categories, the following penalty
#' (Tutz and Gertheiss, 2016),
#' \deqn{J_{\lambda}(\delta)= \sum_{\substack{j=1}}^{p} \sum_{r=1}^{k-2} (\delta_{r+1,j}-\delta_{rj})^{2}}
#' enables the smoothing of response categories such that all
#' category-specific effects associated with the response turn
#' towards a common global effect.
#'
#' @references
#' McCullagh, P (1980). Regression Models for Ordinal Data.
#'     \emph{Journal of the Royal Statistical Society. Series B
#'     (Methodological)}, 42, pp. 109-142.
#'
#' Tutz, G and Gertheiss, J (2016). Regularized Regression
#'     for Categorical Data (With Discussion and Rejoinder).
#'     \emph{Statistical Modelling}, 16, pp. 161-260.
#'
#' @return An object of class "serp" with the components listed below. Other model discription methods include: \code{summary},
#' \code{coef}, \code{predict}, \code{vcov}, \code{anova}, etc.
#' \describe{
#'   \item{aic}{the akaike information criterion.}
#'   \item{bic}{the bayesian information criterion.}
#'   \item{call}{the matched call.}
#'   \item{coef}{a vector of coefficients of the fitted model.}
#'   \item{converged}{a character vector of fit convergence status.}
#'   \item{contrasts}{(where relevant) the contrasts used in the model.}
#'   \item{control}{list of control parameters from \code{serp.control}.}
#'   \item{deviance}{the residual deviance.}
#'   \item{edf}{the (effective) number of degrees of freedom used by the model}
#'   \item{fitted.values}{the fitted probabilities.}
#'   \item{gradient}{a column vector of gradients for the coefficients at the model convergence.}
#'   \item{Hessian}{the hessian matrix for the coefficients at the model convergence.}
#'   \item{iter}{the number of interactions before convergence or non-convergence.}
#'   \item{lambda}{an optimal shrinkage parameter for the penalized slope via minimum residual deviance, cross-validation or user-supplied.}
#'   \item{logLik}{the realized log-likelihood at the model convergence.}
#'   \item{link}{a character vector indicating the link function of the fit.}
#'   \item{model}{model.frame having variables from formula.}
#'   \item{na.action}{(where relevant) information on the treatment of NAs.}
#'   \item{nobs}{the number of observations.}
#'   \item{slope}{a character vector of the type of slope fitted.}
#'   \item{Terms}{the terms structure describing the model.}
#'   \item{tuning}{a character vector specifying the method for choosing an optimal shrinkage parameter.}
#'   \item{value}{the minimum value of the performance metric that yielded the optimal shrinkage parameter (lambda). The brier score is returned for the "cv" tuning, while the residual deviance is reported for the "deviance" tuning.}
#'   \item{ylev}{the number of the response levels.}
#' }
#'
#'
#' @export
#' @examples
#' \dontrun{
#' ## Penalized cumulative logit model of the wine dataset. The optimal
#' ## tuning parameter corresponds to minimum deviance on a continuous
#' ## grid from 0 to a specified maximum value.
#'
#' f1 <- serp(rating ~ temp + contact, slope = "penalize",
#'            reverse = T, link = "logit", tuning = "deviance",
#'            data = wine)
#' summary(f1)
#'
#'
#' ## A user-specified discrete lambda grid for obtainning the optimal
#' ## tuning parameter is allowed for the "cv" and "deviance" tunings.
#'
#' f2 = serp(rating ~ temp + contact, slope = "penalize",
#'           reverse = T, link = "logit", tuning = "deviance",
#'           lambdagrid = seq(0,5,length.out = 50), data = wine)
#' f2
#' head(predict(f2))
#'
#'
#' ## Tuning based on a k-fold cross validation, with k-fold
#' ## ranging from 2 and 10. Default being k = 5.
#'
#' f3 <- serp(rating ~ temp + contact, slope = "penalize",
#'            reverse = F, link = "logit", tuning = "cv",
#'            data = wine)
#' coef(f3)
#'
#'
#' ## The "finite" tuning determines the point along parameter
#' ## shrinkage at which an initially non-finite maximum
#' ## log-likelihood becomes finite. The general model is
#' ## returned by default if log-likelihood initially exists.
#'
#' f4 <- serp(rating ~ temp + contact, slope = "penalize",
#'            reverse = T, link = "logit", tuning = "finite",
#'            data = wine)
#' summary(f4)
#'
#'
#' ## Manual tuning with a user supplied lambda value is also
#' ## possible.
#'
#' f5 <- serp(rating ~ temp + contact, slope = "penalize",
#'            reverse = T, link = "logit", tuning = "manual",
#'            lambda = 43, data = wine)
#' f5
#'
#'
#' ## The unpenalized non-proportional odds model.
#'
#' f6 <- serp(rating ~ temp + contact, slope = "unparallel",
#'            reverse = T, link = "logit", data = wine)
#' summary(f6)
#'
#'
#' ## The unpenalized proportional odds model.
#'
#' f7 <-  serp(rating ~ temp + contact, slope = "parallel",
#'             reverse = T, link = "logit", data = wine)
#' summary(f7)
#'
#'
#' ## The unpenalized partial proportional odds model.
#'
#' f8 <- serp(rating ~ temp + contact, globalvar = ~ temp,
#'            slope = "partial", reverse = T, link = "logit",
#'            data = wine)
#' summary(f8)
#' }
#'
serp <- function(formula,
                 link = c("logit", "probit","loglog", "cloglog", "cauchit"),
                 slope = c("penalize", "parallel", "unparallel", "partial"),
                 tuning = c("deviance","cv","manual", "finite"),
                 reverse = FALSE,
                 lambdagrid = NULL,
                 globalvar=NULL,
                 data,
                 subset=NULL,
                 weights=NULL,
                 weight.type = c("analytic", "frequency"),
                 na.action=NULL,
                 lambda = NULL,
                 contrasts = NULL,
                 control = list(),
                 ...)
{
  function.name <- "serp"
  argnames <- names(list(...))
  mc <- match.call()
  sc <- sys.call()
  checkArg(mc, sc, argnames)
  control <- do.call("serp.control", control)
  m <- match.call(expand.dots = FALSE)
  slope <- match.arg(slope)
  link <- match.arg(link)
  tuning <- match.arg(tuning)
  weight.type <- match.arg(weight.type)
  if(slope=="penalize" && tuning=="manual"){
    if(is.null(lambda))
      stop("user supplied lambda value needed for manual tuning.")
    if(!is.numeric(lambda) | lambda<0 )
      stop("lambda should be numeric and non-negative.")
  }
  if(is.matrix(eval.parent(m$data))) m$data <- as.data.frame(data)
  mf <- match(c("formula", "data", "weights","na.action"), names(m), 0L)
  m <- m[c(1L, mf)]
  m[[1L]] <- quote(stats::model.frame)
  m <- eval(m, parent.frame())
  if (!is.null(subset)) m <- m[subset, , drop=FALSE]
  Terms <- attr(m, "terms")
  y <- as.ordered(model.response(m))
  obs <- length(y)
  if (!is.null(weights)){
    if (!is.numeric(weights))
      stop("weights must be a numeric vector")
    if (any(weights < 0))
      stop("negative weights not allowed")
    if(slope=="penalize" && tuning=="cv" && weight.type!="frequency")
      stop("frequency weights should be used with cv tuning")
    if(weight.type=="frequency"){
      if(all(round(weights) != weights))
        stop("frequency weights must be whole numbers")
      m <- m[rep(seq(nrow(m)), weights), ]
      m <- m[,-ncol(m), drop=FALSE]
      re <- all.vars(Terms)[[1]]
      y <- as.ordered(m[ ,re])
      obs <- length(y)
      wt <- rep(1, obs)
    }else{
      wt <- as.vector(model.weights(m))
    }
  } else wt <- rep(1, obs)
  y <- droplevels(y)
  if(!is.factor(y)) stop("response must be a factor")
  nL <- length(levels(y))
  if(nL <= 2) stop("response must have more than two levels")
  x <- model.matrix(Terms, m, contrasts)
  if(!(is.data.frame(x) | is.matrix(x) | is.numeric(x)))
    stop("x should be a data.frame, matrix or numeric vector")
  if(dim(x)[2]==1){
    vnull <- TRUE
    slope <- "parallel"
    tuning <- "manual"
    lambda <- 0
  } else vnull <- FALSE
  if(is.null(dim(x))){
    x <- cbind(x)
    colnames(x) <- colnames(m)[-1L]
  }
  cons <- attr(x, "contrasts")
  nv <- ifelse(!vnull, dim(x)[2L]-1L, dim(x)[2L])
  if(obs != dim(x)[1L]) stop("variable lengths unequal")
  yMtx <- yMx(y, obs, nL)
  linkf <- lnkfun(link)
  if(slope=='partial'){
    gb <- globalvar <- as.character(all.vars(globalvar))
    if(length(gb)==0) stop("invalid model formula for globalvar")
    if(!length(union(colnames(m), gb)) == length(colnames(m)))
      stop("variable(s) for global effects not recognized")
    if(ncol(x) <= 2L || (ncol(x)-1)==length(gb)) slope <- "parallel"
    if(all.vars(Terms)[[1]] %in% gb)
      stop("model response used as predictor in globalvar")
  }
  xlst <- formxL(x, nL, slope, globalvar, m, vnull)
  yFreq <- colSums(yMtx)/obs
  xMat <- do.call(rbind, xlst)
  if(!vnull) x <- x[ ,-1L, drop=FALSE]
  coln <- colnames(x)
  inibeta <- startv(link, linkf, globalvar, x, yFreq, xMat, nL, nv, slope)
  npar <- length(inibeta)
  nFold <- control$nFold
  useout <- FALSE
  res <- serpfit(lambda, globalvar, x, y, inibeta, xlst, xMat, yMtx, nL,
                 obs, npar, linkf, link, reverse, vnull, control,
                 slope, tuning, nv, nFold, lambdagrid, wt, m)
  coefs <- res[[1]]$coef
  coefs <- est.names(coefs, slope, globalvar,
                     coln, x, m, npar, xMat, nL, vnull, useout)
  ans <- c(list(call=mc, coef= coefs ,model=m, slope=slope,
                globalvar= globalvar,
                tuning=tuning, Terms=Terms, control=control,
                lambda=res[[2]], value=res[[3]]), res[[1]])
  ans$na.action <- attr(m, "na.action")
  ans$contrasts <- cons
  class(ans)<- function.name
  ans
}
