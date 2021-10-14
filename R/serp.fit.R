# SERP fit via Newton-Raphson iteration

serpfit <- function(x, y, wt, yMtx, link, slope, reverse, control,
                    Terms, lambda, lambdaGrid, gridType, tuneMethod,
                    globalEff, cvMetric, mslope, nL, obs,
                    vnull, nvar, m)
{
  if (slope == 'penalize' && !is.null(globalEff))
    slope <- 'partial'
  if (slope == 'partial'){
    if (!is.null(globalEff)) {
      if(class(globalEff) != "formula")
        stop("no object of class formula used in globalEff", call. = FALSE)
      if (grepl('["]', c(globalEff)))
        stop("variable name(s) in quotes not allowed in globalEff",
             call. = FALSE)
      globalEff <- as.character(all.vars(globalEff))
      if (length(globalEff) == 0L)
        stop("wrong input(s) in 'globalEff'", call. = FALSE)
      if (!length(union(colnames(m), globalEff)) == length(colnames(m)))
        stop("unknown variable(s) in globalEff", call. = FALSE)
      if (ncol(x) <= 2L || (ncol(x)-1) == length(globalEff))
        slope <- "parallel"
      if (all.vars(Terms)[[1L]] %in% globalEff)
        stop("response not allowed in 'globalEff'", call. = FALSE)
    } else stop("'globalEff' is unspecified", call. = FALSE)
  }
  link_reversed <- FALSE
  initial_link <- link
  if (link == "cloglog" & reverse == TRUE){
    link_reversed <- TRUE
    link <- "loglog"
    initial_link <- "cloglog"
  }
  if(initial_link == "loglog"){
    if (link == "loglog" & reverse == TRUE){
      link_reversed <- TRUE
      link <- "cloglog"
      initial_link <- "loglog"
    }
  }
  linkf <- lnkfun(link)
  xlst <- formxL(x, nL, slope, globalEff, m, vnull)
  yFreq <- colSums(yMtx)/obs
  xMat <- do.call(rbind, xlst)
  colnames.x <- colnames(x)
  useout <- FALSE
  if (!vnull) x <- x[ ,-1L, drop = FALSE]
  stv <- startv(link, linkf, globalEff, x, yFreq, xMat, nL, nvar, slope)
  npar <- length(stv)
  startval <- est.names(c(stv), slope, globalEff, colnames.x, x, m, npar,
                        xMat, nL, vnull, useout)
  partial.names <- FALSE
  if (mslope == 'penalize' && slope != 'parallel' &&
      !is.null(globalEff)){
    slope <- "penalize"
    partial.names <- TRUE
  }
  if(is.null(lambdaGrid)) lambdaGrid <- control$int.lambdaGrid
  if(!is.numeric(lambdaGrid) || any(is.na(lambdaGrid)) ||
     any(lambdaGrid < 0L) || !(length(lambdaGrid) > 1L))
    stop("lambdaGrid must be a non-negative numeric vector of length > 1",
         call. = FALSE)
  if (any(lambdaGrid > control$lambda.upper.limit))
    stop("values of lambda above 1e10 are not allowed", call. = FALSE)
  if (slope == "penalize"  && tuneMethod == "user"){
    if (is.null(lambda))
      stop("user-supplied lambda value is required for 'user' tuning",
           call. = FALSE)
    if (!is.numeric(lambda) || length(lambda) != 1 || any(lambda < 0))
      stop("lambda should be a single numeric and non-negative value",
           call. = FALSE)
    if (lambda > control$lambda.upper.limit)
      stop("values of lambda above 1e10 are not allowed", call. = FALSE)
  }
  ans <- list(model = m)
  if(slope == "penalize"){
    switch(
      tuneMethod,
      aic = {
        switch(
          gridType,
          discrete = {
            ml <- suppressWarnings(
              vapply(lambdaGrid, dvfun, numeric(1), globalEff, x, y, startval,
                     xlst, xMat, yMtx, nL, obs, npar, linkf, link, vnull,
                     control, slope, wt, tuneMethod, m, mslope, Terms))
            hh <- cbind(lambdaGrid, ml)
            hh <- as.numeric(hh[which.min(hh[,2L]), ])
            minL <- list(minimum = hh[1L], objective = hh[2L])
            lam <- minL$minimum
            ans$lambdaGrid <- lambdaGrid
          },
          fine = {
            minL <- suppressWarnings(
              optimize(dvfun, c(0L,control$maxpen), globalEff,
                       x, y, startval,xlst, xMat, yMtx, nL, obs, npar, linkf,
                       link, vnull,control, slope, wt, tuneMethod, m,
                       mslope, Terms))
            lam <- minL$minimum
          })
      },
      cv = {
        nrFold <- control$nrFold
        switch(
          gridType,
          discrete = {
            ml <- try(suppressWarnings(
              vapply(lambdaGrid, cvfun, numeric(1), x, y, nrFold, linkf, link,
                     m, slope, globalEff, nvar, reverse, vnull, control, wt,
                     cvMetric, mslope, Terms, xlst = xlst, yMtx = yMtx,
                     obs)), silent = TRUE)
            if (inherits(ml, "try-error"))
              stop("cv tuning did not succeed, try using a different tuneMethod")
            hh <- cbind(lambdaGrid, ml)
            hh <- as.numeric(hh[which.min(hh[,2L]), ])
            minL <- list(minimum = hh[1L], objective = hh[2L])
            lam <- minL$minimum
            ans$lambdaGrid <- lambdaGrid
          },
          fine = {
            minL <- try(suppressWarnings(
              optimize(cvfun, c(0L,control$maxpen), x, y, nrFold, linkf,
                       link, m, slope, globalEff, nvar, reverse, vnull,
                       control, wt, cvMetric, mslope, Terms, xlst = xlst,
                       yMtx = yMtx, obs)), silent = TRUE)
            if (inherits(minL, "try-error"))
              stop("cv tuning did not succeed, try using a different tuneMethod")
            lam <- minL$minimum
          })
        ans$nrFold <- nrFold
        ans$cvMetric <- cvMetric
      },
      finite = {
        lx <- dvfun(lambda = 0L, globalEff, x, y, startval,xlst, xMat,
                    yMtx, nL, obs, npar, linkf, link, vnull,control,
                    slope, wt, tuneMethod, m, mslope, Terms)
        if(is.na(lx)){
          minL <- suppressWarnings(
            optimize(f = dvfun, interval = c(0,control$maxpen), globalEff,
                     x, y, startval, xlst, xMat, yMtx, nL, obs, npar, linkf,
                     link, vnull,control, slope, wt, tuneMethod, m,
                     mslope, Terms))
          lam <- minL$minimum
        }else{
          lam <- 0L
          minL <- NULL
        }},
      user = {
        lam <- lambda
        minL <- NULL
        ans$lambda <- lambda
      })
  }else{
    lam <- 0L
    minL <- NULL
  }
  res <- serp.fit(lam, x, y, wt, startval, xlst, xMat,
                  yMtx, nL, obs, npar, linkf, link, vnull,
                  control, slope, globalEff, m, mslope, Terms)
  nmsv <- names(startval)
  misc <- list(colnames.x = colnames.x, colnames.xMat = nmsv, npar = npar,
               variable.null = vnull, convg.no = res$conv,  no.var = nvar)
  delta <- res$coef
  if (reverse) delta <- -1L * delta
  fv <- as.data.frame(res$exact.pr)
  colnames(fv) <- levels(y)
  nfinite <- "non-finite log-likelihood persists, "
  if (!is.null(globalEff)) ans$globalEff <- globalEff
  if (slope == "penalize"){
    if (is.na(res$loglik)){
      if (tuneMethod=="user")
        warning(nfinite,"apply a stronger lambda value or a different ",
                "tuning method")
      if (tuneMethod=="aic" || tuneMethod=="cv"){
        if (!is.null(lambdaGrid))
          warning(nfinite,"increase lambdaGrid upper limit or apply a ",
                  "different tuning method")
        else
          warning(nfinite,"try a different tuning method.")
      }
    }
    if (tuneMethod == "aic" || tuneMethod == "cv"){
      if (tuneMethod == "cv"){
        ans$testError <- if (is.na(res$loglik)) NA else minL$objective
      } else ans$value <- as.numeric(minL$objective)
      if (gridType == "discrete"){
        misc$gridType <- gridType
        misc$grid.range <- if(!is.null(lambdaGrid))
          c(0L, max(lambdaGrid)) else c(0L, control$maxpen)
      } else {
        misc$grid.range <- c(0L, control$maxpen)
        ans$gridType <- gridType
      }
    }
    ans$lambda <- if (is.na(res$loglik)) NA else lam
    if (tuneMethod == "finite") ans$value <- res$loglik
  }
  ans$tuneMethod <- tuneMethod
  if (!vnull){
    hes <- res$info
    gra <- res$score
    dimnames(hes) <- list(names(startval), names(startval))
    dimnames(gra) <- list(names(startval), "gradient")
  } else {
    delta <- delta[nL] + delta[seq_len(nL-1L)]
    hes <- res$info[-nL, -nL]
    gra <- res$score[-nL]
    npar <- nL - 1
    new.name <- names(startval)[seq_len(nL)-1L]
    dimnames(hes) <- list(new.name, new.name)
    names(gra) <- new.name
  }
  if (link_reversed) link <- initial_link
  delta <- as.numeric(delta)
  sl <- if (partial.names) "partial" else slope
  delta <- est.names(delta, slope = sl, globalEff, colnames.x, x, m,
                     npar, xMat, nL, vnull, useout)
  ans$coef <- delta
  ans$logLik <- as.numeric(res$loglik)
  ans$deviance <- -2*ans$logLik
  ans$aic <- ans$deviance + 2*res$edf
  ans$bic <- ans$deviance + log(obs)*res$edf
  ans <- c(list(link = link, ylev = nL, nobs = obs,
                gradient = gra, hessian = hes, fitted.values = fv,
                slope = slope, Terms = Terms, control = control,
                reverse = reverse, converged = res$converged,
                iter = res$iter, message = res$message, misc = misc,
                edf = res$edf, rdf = res$rdf), ans)
  ans
}

serp.fit <- function(lambda, x, y, wt, startval, xlst,
                     xMat, yMtx, nL, obs, npar, linkf, link,
                     vnull, control, slope, globalEff, m, mslope,
                     Terms, xtrace = TRUE)
{
  iter <- 0
  maxits <- control$maxits
  eps <- control$eps
  if ((obs*(nL-1L)) < npar)
    stop("There are ", npar, " parameters but only ", (obs*(nL-1L)),
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
                  slope, m, globalEff, mslope, Terms)
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
      min.ev <- try(min(eigen(info, symmetric = TRUE,
                              only.values = TRUE)$values), silent = TRUE)
      inflation.factor <- 1
      if(inherits(min.ev, "try-error"))
        stop("\nnon-finite eigen values in iteration process", call. = FALSE)
      inflector <- abs(min.ev) + inflation.factor
      info <- info + diag(inflector, nrow(info))
      if(control$trace > 0 && xtrace)
        Trace(iter, maxGrad, obj, delta, step, score, eigen,
              info, trc, inflector, first = (iter == 1), half.iter,
              step.type = "adjust")
      cho <- try(chol(info), silent=TRUE)
      if(inherits(cho, "try-error"))
        stop(gettextf("Cannot compute Newton step at iteration %d",
                      iter), call. = FALSE)
      adj.iter <- adj.iter + 1L
    } else adj.iter <- 0L
    if(adj.iter >= control$maxAdjIter) {
      conv <- 4L
      break
    }
    step <- backsolve(cho, backsolve(cho, score, transpose = TRUE))
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
    if(iter == 2L && fvalues$negprob){
      loglik <- NA
      conv <- 5L
      fvalues <- fvaluesOld
      delta <- deltaOld
      nonfinite <- 1L
      break
    }
    converged <- abs.conv
    half.iter <- 0L
    delta.old <- delta
    obj.old <- obj
    fvalues.old <- fvalues
    while (obj < objOld && (conv == 2L)) {
      delta <- (delta + deltaOld) * 0.5
      fvalues <- prlg(delta, xMat, obs, yMtx, penx, linkf, control, wt)
      pr <- fvalues$pr[,-nL]
      obj <- loglik <- fvalues$logL
      half.iter <- half.iter + 1
      if(control$trace > 0 && xtrace)
        Trace(iter, maxGrad, obj, delta, step, score, eigen,
              info, trc, inflector, first=(iter == 1L), half.iter,
              step.type = "modify")
      if(half.iter >= control$max.half.iter)
      {
        conv <- 3L
        break
      }
    }
    NotImproved <- objOld <= obj
    if (NotImproved) {
      delta <- delta.old
      loglik <- obj.old
      fvalues <- fvalues.old
    }
  }
  if (maxits > 1L && iter >= maxits){
    conv <- 2L
    warning("convergence not obtained in ", maxits,
            " Newton's iterations")
  }
  if (nonfinite && !slope == "penalize")
    warning(control$msg[as.character("s")])
  msg <- control$msg[as.character(conv)][[1L]]
  if(conv <= 1L && control$trace > 0L && xtrace) {
    message("\nSuccessful convergence!\n", msg)
  }
  if(conv > 1 && control$trace > 0L && xtrace) {
    message("\nOptimization failed!\n", msg)
  }
  df <- df.serp(fvalues, xMat, penx, nL)
  res <- c(list(coef = delta, loglik = loglik, info = info, score = score,
                converged = converged, conv = conv, iter = iter,
                message = msg), edf = df$edf, rdf = df$rdf, fvalues)
  return(res)
}
