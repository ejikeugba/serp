# SERP fit via Newton-Raphson iteration
serpfit <- function(lambda, globalvar, x, y, startval, xlst, xMat,
                    yMtx, nL, obs, npar, linkf, link, reverse,
                    vnull,control, slope, tuning, nv, nFold, lambdagrid,
                    wt, m, cverror, mslope, coln, useout, Terms)
{
  if (mslope == 'penalize' && slope != 'parallel' &&
      !is.null(globalvar)) slope <- "penalize"
  if(slope=="penalize"){
    switch(
      tuning,
      deviance = {
        if(!is.null(lambdagrid)){
          ml <- suppressWarnings(
            sapply(lambdagrid, dvfun, globalvar, x, y, startval, xlst, xMat,
                   yMtx, nL, obs, npar, linkf, link, reverse, vnull,
                   control, slope, wt, tuning, m, mslope, Terms))
          hh <- cbind(lambdagrid, ml)
          hh <- hh[which.min(hh[,2L]), ]
          minL <- list(minimum = hh[1L], objective = hh[2L])
          lam <- as.numeric(minL$minimum)
        }else{
          minL <- suppressWarnings(
            optimize(f=dvfun, interval=c(0,control$maxpen), globalvar,
                     x, y, startval,xlst, xMat, yMtx, nL, obs, npar, linkf,
                     link, reverse, vnull,control, slope, wt, tuning, m,
                     mslope, Terms))
          lam <- minL$minimum
        }},
      cv = {
        if(!is.null(lambdagrid)){
          ml <- try(suppressWarnings(
            sapply(lambdagrid, cvfun, x, y, nFold, linkf, link, m, slope,
                   globalvar, nv, reverse, vnull, control, wt, cverror,
                   mslope, tuning, coln, useout, Terms)), silent = TRUE)
          if (inherits(ml, "try-error")) stop("bad input in cv function")
          hh <- cbind(lambdagrid, ml)
          hh <- hh[which.min(hh[,2L]), ]
          minL <- list(minimum = hh[1L], objective = hh[2L])
          lam <- as.numeric(minL$minimum)
        }else{
          minL <- try(suppressWarnings(
            optimize(f=cvfun,interval = c(0,control$maxpen), x, y, nFold,
                     linkf, link, m, slope, globalvar, nv, reverse,
                     vnull, control, wt, cverror, mslope, tuning, coln,
                     useout, Terms)), silent = TRUE)
          if (inherits(minL, "try-error")) stop("bad input in cv function")
          lam <- minL$minimum
        }},
      finite = {
        lx <- dvfun(lambda=0, globalvar, x, y, startval,xlst, xMat,
                    yMtx, nL, obs, npar, linkf, link, reverse,
                    vnull,control, slope, wt, tuning, m, mslope, Terms)
        if(is.na(lx)){
          minL <- suppressWarnings(
            optimize(f=dvfun, interval=c(0,control$maxpen), globalvar,
                     x, y, startval, xlst, xMat, yMtx, nL, obs, npar, linkf,
                     link, reverse, vnull,control, slope, wt, tuning, m,
                     mslope, Terms))
          lam <- minL$minimum
        }else{
          lam <- 0
          minL <- NULL
        }},
      manual = {
        lam <- lambda
        minL <- NULL
      })
  }else{
    lam <- 0
    minL <- NULL
  }
  res <- serp.fit(lam, globalvar, x, y, startval, xlst, xMat,
                  yMtx, nL, obs, npar, linkf, link, reverse,
                  vnull, control, slope, wt, m, mslope, tuning,
                  Terms)
  if (slope == "penalize"){
    if (is.na(res$loglik)){
      if (tuning=="manual")
        warning("non-finite log-likelihood persists, ",
                "try larger values of lambda.")
      if (tuning=="deviance" || tuning=="cv"){
        if (!is.null(lambdagrid))
          warning("non-finite log-likelihood persists, ",
                  "try increasing lambdagrid upper limit.")
        else
          warning("non-finite log-likelihood persists, ",
                  "try other tuning methods.")
      }
    }
    if (tuning=="finite") minL$objective <- res$loglik
  }
  delta <- res$coef
  if (reverse){
    delta <- c(delta[(nL-1):1], delta[nL:nrow(delta)])
    fv <- res$exact.pr[ ,nL:1]
  } else fv <- res$exact.pr
  fv <- as.data.frame(fv)
  colnames(fv) <- levels(y)
  if (!vnull){
    hes <- res$info
    gra <- res$score
  } else {
    delta <- delta[nL] + delta[1:(nL-1)]
    hes <- res$info[-nL, -nL]
    gra <- res$score[-nL]
    npar <- nL - 1
  }
  dv <- c(-2*res$loglik)
  ac <- dv + 2*npar
  bc <- dv + log(obs)*npar
  ans <- list(lambda=lam, value=minL$objective, coef = c(delta),
              logLik = c(res$loglik), deviance = dv, aic = ac, bic = bc,
              converged = res$converged, edf = npar, iter = res$iter,
              ylev = nL, link = link, conv = res$conv, nobs = obs,
              gradient = gra, hessian = hes, fitted.values = fv,
              message = res$message)
  ans
}

serp.fit <- function(lambda, globalvar, x, y, startval, xlst,
                     xMat, yMtx, nL, obs, npar, linkf, link, reverse,
                     vnull, control, slope, wt, m, mslope, tuning,
                     Terms, xtrace=TRUE)
{
  iter <- 0
  maxits <- control$maxits
  eps <- control$eps
  if ((obs*(nL-1)) < npar)
    stop("There are ", npar, " parameters but only ", (obs*(nL-1)),
         " observations")
  conv <- 2L
  delta <- startval
  trc <- control$trace
  converged <- FALSE
  adj.iter <- abs.iter <- nonfinite <- half.iter <- 0L
  while(!converged && iter < maxits)
  {
    iter <- iter+1
    penx <- PenMx(lamv = lambda, delta, nL,
                  slope, m, globalvar, mslope, tuning, Terms)
    fvalues <- prlg(delta, xMat, obs, yMtx, penx, linkf, control, wt)
    pr <- fvalues$pr[,-nL]
    obj <- fvalues$logL
    if(!is.finite(obj))
      stop("Non-finite log-likelihood at starting value")
    SI <- ScoreInfo(x, y, pr, wt, nL, yMtx, xlst, penx, linkf)
    score <- SI$score
    maxGrad <- max(abs(score))
    info <- SI$info
    fvaluesOld <- fvalues
    deltaOld <- delta
    objOld <- obj
    cho <-  try(chol(info), silent = TRUE)
    if(inherits(cho, "try-error")) {
      min.ev <- try(min(eigen(info, symmetric=TRUE,
                              only.values=TRUE)$values), silent=TRUE)
      inflation.factor <- 1
      if(inherits(min.ev, "try-error"))
        stop("\nnon-finite eigen values in iteration path", call.=FALSE)
      inflector <- abs(min.ev) + inflation.factor
      info <- info + diag(inflector, nrow(info))
      if(control$trace > 0 && xtrace)
        Trace(iter, maxGrad, obj, delta, step, score, eigen,
              info, trc, inflector, first=(iter==1), half.iter, step.type = "adjust")
      cho <- try(chol(info), silent=TRUE)
      if(inherits(cho, "try-error"))
        stop(gettextf("Cannot compute Newton step at iteration %d",
                      iter), call.=FALSE)
      adj.iter <- adj.iter + 1L
    } else adj.iter <- 0L
    if(adj.iter >= control$maxAdjIter) {
      conv <- 4L
      break
    }
    step <- backsolve(cho, backsolve(cho, score, transpose=TRUE))
    rel.conv <- (max(abs(step)) < control$relTol)
    delta <- delta + step
    fvalues <- prlg(delta, xMat, obs, yMtx, penx, linkf, control, wt)
    pr <- fvalues$pr[,-nL]
    obj <- loglik <- fvalues$logL
    abs.conv <- eval(control$stopcrit)
    if(abs.conv && !rel.conv) conv <- 1L
    if(abs.conv && rel.conv)  conv <- 0L
    if(control$trace > 0 && xtrace)
      Trace(iter, maxGrad, obj, delta, step, score, eigen, info, trc,
            inflector, first=(iter==1), half.iter, step.type = "full")
    if(iter==2L && fvalues$negprob){
      loglik <- NA
      conv <- 5L
      fvalues <- fvaluesOld
      delta <- deltaOld
      nonfinite <- 1L
      break
    }
    converged <- abs.conv
    half.iter <- 0L
    while (obj < objOld && (conv == 2L)) {
      delta <- (delta + deltaOld) * 0.5
      fvalues <- prlg(delta, xMat, obs, yMtx, penx, linkf, control, wt)
      pr <- fvalues$pr[,-nL]
      obj <- loglik <- fvalues$logL
      half.iter <- half.iter + 1
      if(control$trace > 0 && xtrace)
        Trace(iter, maxGrad, obj, delta, step, score, eigen,
              info, trc, inflector, first=(iter==1), half.iter, step.type = "modify")
      if(half.iter >= control$max.half.iter)
      {
        conv <- 3L
        break
      }
    }
    if(conv == 3L) break
    Improved <- objOld <= obj
    if (!Improved) {
      delta <- deltaOld
      loglik <- objOld
    }
  }
  if (maxits > 1L && iter >= maxits){
    conv <- 2L
    warning("convergence not obtained in ", maxits,
            " Newton-Raphson iterations")
  }
  if (nonfinite && !slope == "penalize")
    warning(control$msg[as.character("s")])
  msg <- control$msg[as.character(conv)][[1L]]
  if(conv <= 1L && control$trace > 0L && xtrace) {
    cat("\n\nSuccessful convergence! ", msg, fill = TRUE)
  }
  if(conv > 1 && control$trace > 0L && xtrace) {
    cat("\n\n Optimization failed!\n", msg, fill = TRUE)
  }
  res <- c(list(coef = delta, loglik = loglik, info = info,
                score = score, converged = converged,
                conv = conv, iter = iter, message = msg), fvalues)
  return(res)
}
