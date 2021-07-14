#' Smooth Effects on Response Penalty for CLM
#'
#' @description Fits cumulative link models (CLMs) with the
#' smooth-effect-on-response penalty (SERP) via a modified Newton-Raphson
#' algorithm. SERP enables the regularization of the parameter space between
#' the general and the restricted cumulative models, with a resultant shrinkage
#' of all subject-specific effects to global effects. The minimum deviance
#' and CV tuning among others, provide the means of arriving at an optimal
#' model in a situation where a user-supplied tuning value is not available.
#' The slope argument allows for the selection of a penalized, unparallel,
#' parallel, or partial slope.
#'
#' @usage serp(
#'      formula,
#'      link = c("logit", "probit","loglog", "cloglog", "cauchit"),
#'      slope = c("penalize", "parallel", "unparallel", "partial"),
#'      tuneMethod = c("deviance", "cv", "finite", "user"),
#'      reverse = FALSE,
#'      lambdaGrid = NULL,
#'      cvMetric = c("brier", "logloss", "misclass"),
#'      gridType = c("discrete", "fine"),
#'      globalEff = NULL,
#'      data,
#'      subset,
#'      weights = NULL,
#'      weight.type = c("analytic", "frequency"),
#'      na.action = NULL,
#'      lambda = NULL,
#'      contrasts = NULL,
#'      control = list(),
#'      ...)
#' @param formula regression formula of the form: response ~ predictors. The
#'   response should be a factor (ordered).
#' @param link sets the link function for the cumulative link model including:
#'   logit, probit, complementary log-log, cloglog, cauchit.
#' @param slope selects the form of coefficients used in the model, with
#'   \code{penalize} denoting the penalized coefficients, \code{unparallel},
#'   \code{parallel} and \code{partial} denoting the unpenalized non-parallel,
#'   parallel and semi-parallel coefficients respectively.
#' @param tuneMethod sets the method of choosing an optimal shrinkage
#'   parameter, including: \code{deviance}, \code{cv}, \code{finite} and
#'   \code{user}. i.e., the lambda value along parameter shrinkage path at
#'   which the fit's residual deviance or the cross-validated test error is
#'   minimal. The finite tuning is used to obtain the model along parameter
#'   shrinkage for which the log-Likelihood exist (is finite). The 'user'
#'   tuning supports a user-supplied lambda value.
#' @param reverse false by default, when true the sign of the linear predictor
#'   is reversed.
#' @param lambdaGrid optional user-supplied lambda grid for the cv and deviance
#'   tuning methods, when the discrete \code{gridType} is chosen. Negative
#'   range of values are not allowed. A short lambda grid could increase
#'   computation time assuming large number of predictors and cases in the
#'   model.
#' @param cvMetric sets the performance metric for the cv tuning, with the
#'   brier score used by default.
#' @param gridType chooses if a discrete or a continuous lambda grid should be
#'   used to select the optimal tuning parameter. The former is used by default
#'   and could be adjusted as desired in \code{serp.control}. The latter
#'   is on the range (0, \code{maxPen}). A user-supplied grid is also possible,
#'   which automatically overrides the internal grid.
#' @param globalEff specifies variable(s) to be assigned global effects during
#'   penalization or when \code{slope} is set to \code{partial}. Variables are
#'   specified as a formula with an empty left hand side, for instance,
#'   globalEff = ~predictors.
#' @param data optional dataframe explaining the variables used in the formula.
#' @param subset specifies which subset of the rows of the data should be used
#'   for fit. All observations are used by default.
#' @param weights optional case weights in fitting. Negative weights are not
#'   allowed. Defaults to 1.
#' @param weight.type chooses between analytic and frequency weights with the
#'   former used by default. The latter should be used when weights are mere
#'   case counts used to compress the data set.
#' @param na.action a function to filter missing data.
#' @param lambda a user-supplied single numeric value for the tuning parameter
#'   when using the \code{user} tuning method. Negative values are not
#'   allowed.
#' @param contrasts a list of contrasts to be used for some or all of the
#'   factors appearing as variables in the model formula.
#' @param control A list of fit control parameters to replace default values
#'   returned by \code{serp.control}. Values not set assume default values.
#' @param ... additional arguments.
#' @importFrom ordinal pgumbel
#' @importFrom ordinal dgumbel
#' @importFrom ordinal qgumbel
#' @importFrom stats BIC
#' @importFrom stats coef
#' @importFrom stats dcauchy
#' @importFrom stats dlogis
#' @importFrom stats dnorm
#' @importFrom stats formula
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats model.weights
#' @importFrom stats na.omit
#' @importFrom stats optimize
#' @importFrom stats pcauchy
#' @importFrom stats pchisq
#' @importFrom stats plogis
#' @importFrom stats pnorm
#' @importFrom stats printCoefmat
#' @importFrom stats qcauchy
#' @importFrom stats qlogis
#' @importFrom stats qnorm
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
#' of regression parameters and \eqn{F} a continuous distribution
#' function. This model assumes that the effect of \eqn{x} does not
#' depend on the  category. However, with this assumption relaxed,
#' one obtains the following general cumulative model:
#'
#' \deqn{P(Y\leq r|x) = F(\delta_{0r} + x^T\delta_{r}),}
#'
#' where r=1,\dots,k-1. This model, however, has the stochastic ordering
#' property, which implies that \eqn{P(Y\leq r-1|x) < P(Y\leq r|x)}
#' holds for all \eqn{x} and all categories \eqn{r}. Such assumption
#' is often problematic, resulting in unstable likelihoods with
#' ill-conditioned parameter space during the iterative procedure.
#'
#' SERP offers a means of arriving at stable estimates of the general model.
#' It provides a form of regularization that is based on minimizing the
#' penalized log-likelihood:
#'
#' \deqn{l_{p}(\delta)=l(\delta)-J_{\lambda}(\delta)}
#'
#' where \eqn{l(\delta)}, is the log-likelihood of the general cumulative
#' model and \eqn{J_{\lambda}(\delta)=\lambda J(\delta)} the penalty
#' function weighted by the turning parameter \eqn{\lambda}. Assuming an
#' ordered categorical outcome \eqn{Y \in \{1,\dots,k\}}, and considering
#' that the corresponding parameters \eqn{\delta_{1j},\dots \delta_{k-1,j}}
#' vary smoothly over the categories, the following penalty
#' (Tutz and Gertheiss, 2016),
#'
#' \deqn{J_{\lambda}(\delta)= \sum_{j=1}^{p} \sum_{r=1}^{k-2}
#' (\delta_{r+1,j}-\delta_{rj})^{2}}
#'
#' enables the smoothing of response categories such that all
#' category-specific effects associated with the response turn towards a
#' common global effect. SERP could also be applied to a semi-parallel model
#' with only the category-specific part of the model penalized.
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
#' @return An object of class \code{serp} with the components listed below,
#' depending on the type of slope modeled. Other summary methods include:
#'  \code{summary}, \code{coef}, \code{predict}, \code{vcov},
#' \code{anova}, \code{errorMetrics}, etc.
#'
#' \describe{
#'   \item{aic}{the akaike information criterion.}
#'   \item{bic}{the bayesian information criterion.}
#'   \item{call}{the matched call.}
#'   \item{coef}{a vector of coefficients of the fitted model.}
#'   \item{converged}{a character vector of fit convergence status.}
#'   \item{contrasts}{(where relevant) the contrasts used in the model.}
#'   \item{control}{list of control parameters from \code{serp.control}.}
#'   \item{cvMetric}{the performance metric used for cv tuning.}
#'   \item{deviance}{the residual deviance.}
#'   \item{edf}{the (effective) number of degrees of freedom used by the model}
#'   \item{fitted.values}{the fitted probabilities.}
#'   \item{globalEff}{variable(s) in model treated as global effect(s)}
#'   \item{gradient}{a column vector of gradients for the coefficients at the
#'         model convergence.}
#'   \item{Hessian}{the hessian matrix for the coefficients at the model
#'         convergence.}
#'   \item{iter}{number of interactions before convergence or non-convergence.}
#'   \item{lambda}{a user-supplied single numeric value for the \code{user}
#'         tuning tuning method.}
#'   \item{lambdaGrid}{a numeric vector of lambda values used to determine the
#'         optimum tuning parameter.}
#'   \item{logLik}{the realized log-likelihood at the model convergence.}
#'   \item{link}{character vector indicating the link function of the fit.}
#'   \item{message}{character vector stating the type of convergence obtained}
#'   \item{misc}{a list to hold miscellaneous fit information.}
#'   \item{model}{model.frame having variables from formula.}
#'   \item{na.action}{(where relevant) information on the treatment of NAs.}
#'   \item{nobs}{the number of observations.}
#'   \item{nrFold}{the number of k-fold cross validation for the cv tuning
#'         method. Default to k = 5.}
#'   \item{reverse}{a logical vector indicating the the direction of the
#'         cumulative probabilities. Default to P(Y<=r).}
#'   \item{slope}{a character vector indicating the type of slope parameters
#'         fitted. Default to \code{penalize}.}
#'   \item{Terms}{the terms structure describing the model.}
#'   \item{testError}{numeric value of the cross-validated test error at which
#'         the optimal tuning parameter emerged.}
#'   \item{tuneMethod}{a character vector specifying the method for choosing an
#'         optimal shrinkage parameter.}
#'   \item{value}{numeric value of the deviance or the minus log-likelihood of
#'         the optimal model for the 'deviance' and the 'finite' tuning methods
#'         respectively.}
#'   \item{ylev}{the number of the response levels.}
#' }
#' @export
#' @examples
#'
#' ## The unpenalized non-proportional odds model. (with Unbounded estimates)
#' f1 <- serp(rating ~ temp + contact, slope = "unparallel",
#'            reverse = TRUE, link = "logit", data = wine)
#' coef(f1)
#' logLik(f1)
#'
#' ## The penalized non-proportional odds model (with Improved estimates)
#' f2 <- serp(rating ~ temp + contact, slope = "penalize",
#'            link = "logit", reverse = TRUE, tuneMethod = "deviance",
#'            lambdaGrid = 10^seq(1, -1, length.out=5), data = wine)
#' coef(f2)
#' predict(f2, type = "class")
#'
#' ## The unpenalized proportional odds model (with constrained estimates).
#' f3 <-  serp(rating ~ temp + contact, slope = "parallel",
#'             reverse = FALSE, link = "logit", data = wine)
#' summary(f3)
#' confint(f3)
#' errorMetrics(f3)
#'
serp <- function(
  formula,
  link = c("logit", "probit", "loglog", "cloglog", "cauchit"),
  slope = c("penalize", "parallel", "unparallel", "partial"),
  tuneMethod = c("deviance","cv", "finite", "user"),
  reverse = FALSE,
  lambdaGrid = NULL,
  cvMetric = c("brier", "logloss", "misclass"),
  gridType = c("discrete", "fine"),
  globalEff= NULL,
  data,
  subset = NULL,
  weights = NULL,
  weight.type = c("analytic", "frequency"),
  na.action = NULL,
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
  tuneMethod <- match.arg(tuneMethod)
  weight.type <- match.arg(weight.type)
  cvMetric <- match.arg(cvMetric)
  gridType <- match.arg(gridType)
  if (is.matrix(eval.parent(m$data))) m$data <- as.data.frame(data)
  mf <- match(c("formula", "data", "weights","na.action"), names(m), 0L)
  m <- m[c(1L, mf)]
  m[[1L]] <- quote(stats::model.frame)
  m <- eval(m, parent.frame())
  if (!is.null(subset)) {
    if (!is.numeric(subset) || any(subset < 0L) || any(is.na(subset)) ||
        !all(subset == floor(subset)))
      stop("subset indices must be positive whole numbers")
    else m <- m[subset, , drop = FALSE]
  }
  Terms <- attr(m, "terms")
  y <- model.response(m)
  if (is.null(y)) stop("response missing in formula")
  if (!is.ordered(y)) stop("response must be an ordered factor")
  if (length(grep(all.vars(Terms)[1L], mc$formula)) > 1L)
    stop("response not allowed as predictor")
  obs <- length(y)
  if (!is.null(weights)){
    if (!is.numeric(weights) || any(is.na(weights)))
      stop("weights should be numeric vector with no NA's")
    if (any(weights < 0))
      stop("negative weights not allowed")
    if (slope == "penalize" && tuneMethod == 'cv' &&
        weight.type != "frequency")
      stop("frequency weights only is allowed for 'cv' tuning")
    if (weight.type == "frequency"){
      if(any(round(weights) != weights))
        stop("frequency weights must be whole numbers")
      m <- m[rep(seq(nrow(m)), weights), ]
      m <- m[,-ncol(m), drop = FALSE]
      re <- all.vars(Terms)[[1L]]
      y <- as.ordered(m[ ,re])
      obs <- length(y)
      wt <- rep(1L, obs)
    }else{
      wt <- as.vector(model.weights(m))
    }
  } else wt <- rep(1L, obs)
  y <- droplevels(y)
  nL <- length(levels(y))
  if (nL <= 2L) stop("response must have 3 or more levels")
  x <- model.matrix(Terms, m, contrasts)
  if (!(is.data.frame(x) || is.matrix(x) || is.numeric(x)))
    stop("x must be a data.frame, matrix or numeric vector")
  if (dim(x)[2L] == 1L){
    vnull <- TRUE
    slope <- "parallel"
    tuneMethod <- "user"
    lambda <- 0L
  } else vnull <- FALSE
  if (is.null(dim(x))){
    x <- cbind(x)
    colnames(x) <- colnames(m)[-1L]
  }
  cons <- attr(x, "contrasts")
  mslope <- slope
  nvar <- ifelse(!vnull, dim(x)[2L] - 1L, dim(x)[2L])
  if (obs != dim(x)[1L]) stop("variable lengths unequal")
  yMtx <- yMx(y, obs, nL)
  ans <- serpfit(x, y, wt, yMtx, link, slope, reverse, control,
                 Terms, lambda, lambdaGrid, gridType, tuneMethod,
                 globalEff, cvMetric, mslope, nL, obs,
                 vnull, nvar, m)
  ans <- c(list(call = mc), ans)
  ans$na.action <- attr(m, "na.action")
  ans$contrasts <- cons
  class(ans) <- function.name
  ans
}
