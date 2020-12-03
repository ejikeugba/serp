# The serp fit
serpfit <- function(lambda, globalvar, x, y, inibeta, xlst, xMat,
                    yMtx, nL, obs, npar, linkf, link, reverse,
                    vnull,control, slope, tuning, nv, nFold, lambdagrid, wt, m)
{
  if(slope=="penalize"){
    switch(
      tuning,
      deviance = {
        if(!is.null(lambdagrid)){
          ml <- suppressWarnings(
            sapply(lambdagrid, dvfun, globalvar, x, y, inibeta, xlst, xMat,
                   yMtx, nL, obs, npar, linkf, link, reverse, vnull,
                   control, slope, wt, tuning))
          hh <- cbind(lambdagrid, ml)
          hh <- hh[which.min(hh[,2L]), ]
          minL <- list(minimum = hh[1L], objective = hh[2L])
          lam <- as.numeric(minL$minimum)
        }else{
          minL <- suppressWarnings(
            optimize(f=dvfun, interval=c(0,control$maxpen), globalvar,
                     x, y, inibeta,xlst, xMat, yMtx, nL, obs, npar, linkf,
                     link, reverse, vnull,control, slope, wt, tuning))
          lam <- minL$minimum
        }},
      cv = {
        if(!is.null(lambdagrid)){
          ml <- try(suppressWarnings(
            sapply(lambdagrid, cvfun, x, y, nFold, linkf, link, m, slope,
                   globalvar, nv, reverse, vnull, control, wt)), silent = TRUE)
          if (inherits(ml, "try-error")) stop("inherited object in cvfun undefined")
          hh <- cbind(lambdagrid, ml)
          hh <- hh[which.min(hh[,2L]), ]
          minL <- list(minimum = hh[1L], objective = hh[2L])
          lam <- as.numeric(minL$minimum)
        }else{
          minL <- try(suppressWarnings(
            optimize(f=cvfun,interval = c(0,control$maxpen), x, y, nFold,
                     linkf, link, m, slope, globalvar, nv, reverse,
                     vnull, control, wt)), silent = TRUE)
          if (inherits(minL, "try-error")) stop("inherited object in cvfun undefined")
          lam <- minL$minimum
        }},
      finite = {
        lx <- dvfun(lambda=0, globalvar, x, y, inibeta,xlst, xMat,
                    yMtx, nL, obs, npar, linkf, link, reverse,
                    vnull,control, slope, wt, tuning)
        if(is.na(lx)){
          minL <- suppressWarnings(
            optimize(f=dvfun, interval=c(0,control$maxpen), globalvar,
                     x, y, inibeta, xlst, xMat, yMtx, nL, obs, npar, linkf,
                     link, reverse, vnull,control, slope, wt, tuning))
          lam <- minL$minimum
        }else{
          lam <- 0
          minL <- NULL
        }
      },
      manual = {
        lam <- lambda
        minL <- NULL
      }
    )
  }else{
    lam <- 0
    minL <- NULL
  }
  res <- serp.fit(lam, globalvar, x, y, inibeta, xlst, xMat,
                  yMtx, nL, obs, npar, linkf, link, reverse,
                  vnull, control, slope, wt)
  if (slope == "penalize"){
    if (is.na(res$logLik)){
      if (tuning=="manual")
        warning("non-existent likelihood, try larger values of lambda.")
      if (tuning=="deviance" || tuning=="cv"){
        if (!is.null(lambdagrid))
          warning("non-existent likelihood, try increasing lambdagrid upper limit.")
        else
          warning("non-existent likelihood, try other tuning methods.")
      }
    }
    if (tuning=="finite") minL$objective <- res$logLik
  }

  list(res, lambda=lam, objective=minL$objective)
}



serp.fit <- function(lambda, globalvar, x, y, inibeta, xlst,
                     xMat, yMtx, nL, obs, npar, linkf, link, reverse,
                     vnull, control, slope, wt)
{
  iter <- 0
  converged <- FALSE
  while( (!converged) & (iter < control$maxits) )
  {
    iter <- iter+1
    delta <- inibeta
    if((obs*(nL-1)) < npar)
      stop("There are ", npar, " parameters but only ", (obs*(nL-1)), " observations ")
    penx <- PenMx(lamv = lambda, delta, nL, slope)
    problik <- prlg(delta, xMat, obs, yMtx, penx, linkf, control, wt)
    pr <- problik$pr[,-nL]
    loglik <- objOld <- problik$llik
    pf <- cbind(pr, 1-rowSums(pr))
    pf <- pf/rowSums(pf)
    lp <- pf[, nL]
    pr <- pf[, -nL, drop = FALSE]
    g  <- yMtx[, nL]/lp
    g[yMtx[, nL] == 0] <- 0
    yp <- yMtx[, -nL, drop = FALSE]/pr
    yp[yMtx[, -nL] == 0] <- 0
    etaMat <- suppressWarnings(t(apply(t(apply(pr,1,cumsum)), 1, linkf$qfun)))
    uMat <- yp - g
    uMat <- split(uMat, 1:nrow(etaMat))
    if(is.null(attr(yMtx, "wts")))
      wts <- rowSums(yMtx) else wts <- attr(yMtx, "wts")
    pw <- wts/pr
    ff <- wts/lp
    q <- lapply(split(etaMat, 1:nrow(etaMat)), etfun, nL, linkf)
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
    invinfo  <-  try(chol2inv(chol(info)), silent = TRUE)
    if (!inherits(invinfo, "try-error")){
      oldscore <- score
      oldinfo <- info
      oldinvinfo <- invinfo
    } else {
      info <- if (exists("oldinfo")) oldinfo
      else
        stop("Hessian is not positive definite")
      score <- oldscore
      invinfo <- oldinvinfo
      break
    }
    ndelta <- delta + invinfo %*% score
    problik <- prlg(delta=ndelta, xMat, obs, yMtx, penx, linkf, control, wt)
    pr  <- problik$pr[,-dim(problik$pr)[2L]]
    obj <- loglik <- problik$llik
    if (is.nan(obj)) obj <- loglik <- Inf
    inibeta <- ndelta
    dif <- obj - objOld
    converged <- dif < control$eps
  }
  if (!inherits(invinfo, "try-error"))
    newdelta <- ndelta
  else
    stop("model not of full-rank")
  if (control$maxits > 1 && iter >= control$maxits)
    warning("convergence not obtained in ", control$maxits, " iterations")
  if (is.nan(loglik) || problik$negprob){
    loglik <- NA
    if(!slope=="penalize")
      warning("stochastic order assumption resulted in intersecting linear/additive \npredictors. Try using penalized or parallel slopes, or other link functions.")
  }
  if(reverse)
    newdelta <- -1 * c(newdelta)
  else
    newdelta <- c(newdelta)
  if(!vnull){
    hes <- info
    gra <- score
  }else {
    newdelta <- newdelta[nL] + newdelta[1:(nL-1)]
    hes <- info[-nL, -nL]
    gra <- score[-nL]
    npar <- nL - 1
  }
  list(coef= newdelta,
       logLik = c(loglik),
       deviance = -2*problik$llik,
       aic = c(-2*loglik + 2*npar),
       bic = c(-2*loglik + log(obs)*npar),
       gradient = gra,
       hessian = hes,
       problik = problik,
       converged = c(converged),
       edf = npar,
       iter = c(iter),
       fitted.values = problik$pr,
       ylev = nL,
       link = link,
       nobs = obs
  )
}
