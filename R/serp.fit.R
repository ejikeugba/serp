# SERP fit via Newton-Raphson iteration
serpfit <- function(lambda, globalvar, x, y, startval, xlst, xMat,
                    yMtx, nL, obs, npar, linkf, link, reverse,
                    vnull,control, slope, tuning, nv, nFold, lambdagrid,
                    wt, m, cverror, mslope, coln, useout)
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
                   control, slope, wt, tuning, m, mslope))
          hh <- cbind(lambdagrid, ml)
          hh <- hh[which.min(hh[,2L]), ]
          minL <- list(minimum = hh[1L], objective = hh[2L])
          lam <- as.numeric(minL$minimum)
        }else{
          minL <- suppressWarnings(
            optimize(f=dvfun, interval=c(0,control$maxpen), globalvar,
                     x, y, startval,xlst, xMat, yMtx, nL, obs, npar, linkf,
                     link, reverse, vnull,control, slope, wt, tuning, m,
                     mslope))
          lam <- minL$minimum
        }},
      cv = {
        if(!is.null(lambdagrid)){
          ml <- try(suppressWarnings(
            sapply(lambdagrid, cvfun, x, y, nFold, linkf, link, m, slope,
                   globalvar, nv, reverse, vnull, control, wt, cverror,
                   mslope, tuning, coln, useout)), silent = TRUE)
          if (inherits(ml, "try-error")) stop("bad input in cv function")
          hh <- cbind(lambdagrid, ml)
          hh <- hh[which.min(hh[,2L]), ]
          minL <- list(minimum = hh[1L], objective = hh[2L])
          lam <- as.numeric(minL$minimum)
        }else{
          minL <- try(suppressWarnings(
            optimize(f=cvfun,interval = c(0,control$maxpen), x, y, nFold,
                     linkf, link, m, slope, globalvar, nv, reverse,
                     vnull, control, wt, cverror, mslope, tuning, coln, useout)), silent = TRUE)
          if (inherits(minL, "try-error")) stop("bad input in cv function")
          lam <- minL$minimum
        }},
      finite = {
        lx <- dvfun(lambda=0, globalvar, x, y, startval,xlst, xMat,
                    yMtx, nL, obs, npar, linkf, link, reverse,
                    vnull,control, slope, wt, tuning, m, mslope)
        if(is.na(lx)){
          minL <- suppressWarnings(
            optimize(f=dvfun, interval=c(0,control$maxpen), globalvar,
                     x, y, startval, xlst, xMat, yMtx, nL, obs, npar, linkf,
                     link, reverse, vnull,control, slope, wt, tuning, m, mslope))
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
                  vnull, control, slope, wt, m, mslope, tuning)
  if (slope == "penalize"){
    if (is.na(res$logLik)){
      if (tuning=="manual")
        warning("Infinite log-likelihood persists, ",
                "try larger values of lambda.")
      if (tuning=="deviance" || tuning=="cv"){
        if (!is.null(lambdagrid))
          warning("Infinite log-likelihood persists, ",
                  "try increasing lambdagrid upper limit.")
        else
          warning("Infinite log-likelihood persists, ",
                  "try other tuning methods.")
      }
    }
    if (tuning=="finite") minL$objective <- res$logLik
  }
  list(res, lambda=lam, objective=minL$objective)
}

serp.fit <- function(lambda, globalvar, x, y, startval, xlst,
                     xMat, yMtx, nL, obs, npar, linkf, link, reverse,
                     vnull, control, slope, wt, m, mslope, tuning)
{
  iter <- 0
  maxits <- control$maxits
  eps <- control$eps
  converged <- FALSE
  while( (!converged) & (iter < maxits) )
  {
    iter <- iter+1
    delta <- startval
    if ((obs*(nL-1)) < npar)
      stop("There are ", npar, " parameters but only ", (obs*(nL-1)),
           " observations")
    penx <- PenMx(lamv = lambda, delta, nL,
                  slope, m, globalvar, mslope, tuning)
    pfit1 <- prlg(delta, xMat, obs, yMtx, penx, linkf,
                  control, wt)
    pr <- pfit1$pr[,-nL]
    objOld <- pfit1$llik
    um <- uMatFun(pr, yMtx, linkf, nL)
    uMat <- um$uMat
    if (is.null(attr(yMtx, "wts")))
      wts <- rowSums(yMtx) else wts <- attr(yMtx, "wts")
    pw <- wts/um$pr
    ff <- wts/um$lp
    q <- lapply(split(um$etaMat, 1:nrow(um$etaMat)), etfun, nL, linkf)
    invmat <- lapply(split(pw, 1:nrow(pw)), diag)
    Infx <- mapply(function(invt, ff, x, q, wt){
      inv <- invt + ff
      w <- crossprod(q,inv) %*% q
      fm <- wt * crossprod(x, w %*% x)
    }, invmat, ff, xlst, q, wt, SIMPLIFY = FALSE)
    info <- Reduce("+", Infx) + penx$infopen
    sc <- mapply(function(q, uMat, x, wt) wt*crossprod(x, crossprod(q, uMat)),
                 q, uMat, xlst, wt, SIMPLIFY = FALSE)
    score <- Reduce("+", sc) - penx$scorepen
    cholisk  <-  try(chol(info), silent = TRUE)
    if (inherits(cholisk, "try-error"))
      stop("Numerical problems in iteration process")
    newdelta <- delta + backsolve(cholisk,
                                  backsolve(cholisk, score, transpose=TRUE))
    pfit1 <- prlg(delta=newdelta, xMat, obs, yMtx, penx, linkf, control, wt)
    objNew <- pfit1$llik
    loglik <- objNew
    startval <- newdelta
    pfit2 <- prlg(delta=newdelta, xMat, obs, yMtx, penx, linkf,
                  control, wt, intLogL = FALSE)
    converged <- eval(control$stopcrit)
  }
  if (maxits > 1 && iter >= maxits)
    warning("convergence not obtained in ", maxits, " iterations")
  if (is.nan(loglik) || pfit1$negprob){
    loglik <- NA
    if(!slope == "penalize")
      warning("Stochastic ordering assumption unmet. ",
              "Consider using penalized, parallel or ",
              "partial slope, or other link function.")
  }
  if (reverse){
    newdelta <- c(newdelta[(nL-1):1], newdelta[nL:nrow(newdelta)])
    fv <- pfit1$prob[ ,nL:1]
  } else fv <- pfit1$prob
  fv <- as.data.frame(fv)
  colnames(fv) <- levels(droplevels(y))
  if (!vnull){
    hes <- info
    gra <- score
  } else {
    newdelta <- newdelta[nL] + newdelta[1:(nL-1)]
    hes <- info[-nL, -nL]
    gra <- score[-nL]
    npar <- nL - 1
  }
  list(coef= c(newdelta),
       logLik = c(loglik),
       deviance = -2*pfit1$llik,
       aic = c(-2*loglik + 2*npar),
       bic = c(-2*loglik + log(obs)*npar),
       gradient = gra,
       hessian = hes,
       problik = pfit2,
       converged = c(converged),
       edf = npar,
       iter = c(iter),
       fitted.values = fv,
       ylev = nL,
       link = link,
       nobs = obs)
}
