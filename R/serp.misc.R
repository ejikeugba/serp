# Internal functions
formxL <- function(x, nL, slope, globalvar, m, vnull, ...)
{
  Xdt <- function(x, nL)
  {
    nv  <- length(x)
    zs <- rep(0, nv)
    bc <- do.call(rbind, lapply(1:(nL-1), function(i)
      c(rep(zs, i-1), x, rep(zs, nL-1-i))
    ))
    bc
  }
  h1 <- diag(nL-1)
  nv <- dim(x)[2]
  lq <- as.vector(t(matrix(1:(nv * (nL-1)), nv, nL-1)))
  xL <- lapply(1:nrow(x), function(u)
  {
    if (slope == 'partial'){
      gb <- unlist(globalvar)
      subg <- subset(m, select = c(gb))
      subm <- model.matrix(~., data = subg)
      subn <- colnames(subm)[-1L]
      xnames <- colnames(x)
      if (any(grepl(":", xnames))){
        subn2 <- xnames[which(grepl(":", xnames))]
        subn <- c(subn, subn2)
      }
      col.num <- seq_len(ncol(x))
      col.glo <- c(1L, which(xnames %in% subn))
      specifc <- setdiff(col.num, col.glo[-1L])
      if(!dim(x)[1L] == 1L)
        xs <- x[,specifc]
      else
        xs <- subset(x, select = specifc)
      nv <- dim(xs)[2]
      lq <- as.vector(t(matrix(1:(nv * (nL-1)), nv, nL-1)))
      xx <- xs[u, ]
      xL1 <- Xdt(xx, nL)[ ,lq]
      xL1 <- xL1[,-c(1L:nL-1L)]
      cname <- unique(colnames(xL1)[colnames(xL1) != ""])
      colnames(xL1) <- rep(cname, each = nL-1)
      if (!dim(x)[1] == 1L)
        xg <- x[,col.glo]
      else
        xg <- subset(x, select = col.glo)
      nv <- dim(xg)[2]
      xx <- xg[u, ]
      h2 <- rbind(xx)[rep(1, nL-1), ,drop = FALSE]
      if (!vnull)
        xL2 <- h2[,-1L, drop = FALSE]
      else
        xL2 <- h2
      XL3 <- cbind(xL1, xL2)
      xListi <- cbind(h1, XL3[ , order(colnames(XL3))])
    }else{
      xx <- x[u, ]
      if (slope == 'parallel'){
        h2 <- rbind(xx)[rep(1, nL-1), ,drop = FALSE]
        if (!vnull)  xListi <- cbind(h1, h2[,-1L])
        else xListi <- cbind(h1, h2)
      }
      if (slope == 'unparallel'|| slope == 'penalize'){
        xListi <- Xdt(xx, nL)[ ,lq]
      }
      rownames(xListi) <- NULL
      xListi
    }
  })
  xL
}

PenMx <- function(lamv, delta, nL, slope)
{
  nvar <- length(delta)-(nL-1)
  diagx <- rep(-1, nvar-1)
  mt <- diag(-1, nrow = nvar-1, ncol = nvar)
  mt <- 1*((row(mt) == col(mt) - (1)) + 0)
  diag(mt) <- diagx
  r <- 0
  while(r < nvar-1)
  {
    r <- r + nL-1
    if (r > nvar-1) break
    mt[r, ] <- 0
  }
  Dx <- t(mt) %*% mt
  pM <- rbind(matrix(0, nL-1, nvar+nL-1),
              cbind(matrix(0, nvar, nL-1), Dx))
  if (slope=="parallel" || slope=="partial") lamv <- 0
  py1 <- lamv * pM
  py2 <- lamv * pM %*% delta
  py3 <- lamv * t(delta) %*%  pM  %*% delta
  list(infopen=py1, scorepen=py2, logpen=py3)
}

prlg <- function(delta, xMat, obs,
                 yMtx, penx, linkf, control, wt){
  hh <- xMat %*% delta
  invlink <- 1 - apply(hh, 1, linkf$pfun)
  ieta <- cbind(1,matrix(invlink, nrow=obs, byrow=TRUE))
  if(dim(ieta)[1L]==1L)
    pr <- prob <- c(ieta[, -ncol(ieta)] - ieta[, -1], ieta[,dim(ieta)[2L]])
  else
    pr <- prob <- cbind(ieta[, -ncol(ieta)] - ieta[, -1], ieta[,dim(ieta)[2L]] )
  np <- any(pr<0)

  double.eps <- sqrt(.Machine$double.xmin)
  pr[pr <= double.eps] <- double.eps
  llik <- suppressWarnings(sum(wt * yMtx * log(pr)))
  llik2 <- suppressWarnings(sum(wt * yMtx * log(prob)))

  lk <- llik + penx$logpen
  list(pr=pr, llik=lk, negprob=np, prob=prob, llik2=llik2)
}

checkArg <- function (mcall, scall, argnames)
{
  nm <- names(as.list(mcall))
  if (!"formula" %in% nm)
    stop("Model needs a formula", call. = FALSE)
  nd <- list(m1=names(mcall), m2=names(scall))
  nd <- setdiff(nd$m1, nd$m2)
  if (length(nd) > 1) stop("unassigned value in serp function.", call. = FALSE)
  check <- argnames %in% names(formals(cat))
  if (any(!check)) anyerr <- TRUE else anyerr <- FALSE
  if (length(check)>1)
    err <- sprintf("unused arguments: \"%s\"",
                   paste(argnames[!check], collapse = ", "))
  else
    err <- sprintf("unused argument: \"%s\"",
                   paste(argnames[!check], collapse = ", "))
  if (anyerr) stop(err, call. = FALSE)
}

startv <- function(link, linkf, globalvar, x, yFreq, xMat, nL, nv, slope)
{
  ut <- linkf$qfun(cumsum(yFreq[-nL]))
  if (link == "cauchit") ut <- tan(pi*(cumsum(yFreq[-nL]) - 0.5))
  if (slope == "parallel") ibeta <- c(ut, rep(0, nv))
  if (slope == "partial") ibeta <- c(ut, rep(0, ncol(xMat)-(nL-1L)))
  if (slope == "unparallel" || slope == "penalize") ibeta <- c(ut, rep(0, nv * (nL - 1)))
  ibeta
}

est.names <- function(coef, slope, globalvar,
                      coln, x, m, npar, xMat, nL, vnull, useout){
  if (slope == "parallel"){
    rr <- c(sprintf("(Intercept):%d", 1:(nL-1)), coln)
    if(vnull == TRUE) rr <- rr[1:(nL-1)]
    names(coef) <- rr
  }
  if (slope == "unparallel" || slope == "penalize"){
    xname <- rep(coln, each = nL-1)
    rr <- c(sprintf("(Intercept):%d", 1:(nL-1)),
            paste(xname, 1:(nL-1), sep = ":"))
    if (vnull == TRUE) rr <- rr[1:(nL-1)]
    names(coef) <- rr
  }
  if (slope == "partial"){
    xnam <- colnames(xMat)[-c(1:nL-1)]
    gb <- unlist(globalvar)
    subg <- subset(m, select = c(gb))
    subm <- model.matrix(~., data = subg)
    subn <- colnames(subm)[-1L]
    if (any(grepl(":", xnam))){
      subn2 <- xnam[which(grepl(":", xnam))]
      subn <- c(subn, subn2)
    }
    col.num <- seq_len(length(xnam))
    col.glo <- which(xnam %in% subn)
    specifc <- setdiff(col.num, col.glo)
    ss <- xnam[specifc]
    if (vnull == TRUE) ss <- ss[1:(nL-1)]
    ww <- rep(coln, each=nL-1)
    ww <- ww[!duplicated(ww, incomparables = ss)]
    origvec <- varnames(xnam, coef, nL)
    posvec <- varnames(ww, coef, nL)
    coef <- origvec[names(posvec)]
    if(useout) coef <- posvec[names(origvec)]
  }
  coef
}

varnames <- function(cofnames, coef, nL){
  thresh <- sprintf("(Intercept):%d", 1:(nL-1))
  cfn <- make.unique(as.character(cofnames), sep=":" )
  ind <- !duplicated(cofnames)
  for (i in 1:length(cfn[ind])){
    cfn[ind][i]<-ifelse(cofnames[ind][i] %in% cofnames[duplicated(cofnames)], gsub("(.*)","\\1:0", cfn[ind][i]),
                        cfn[ind][i])
  }
  mx <- sub(".*:([0-9]+)","\\1",grep(":([0-9]*)$",cfn,value=T) )
  mn <- sub("(.*:)[0-9]+","\\1",grep(":([0-9]*)$",cfn,value=T) )
  mxm <- as.numeric(mx)+1
  cfn[grep(":([0-9]*)$",cfn)] <- paste0(mn,mxm)
  names(coef) <- c(thresh, cfn)
  coef
}

levformat <- function (pr, digits)
  paste(format(pr * 1e2, scientific = FALSE, trim = F, digits = digits),"%")

etfun <- function(eta, nL, linkf)
{
  mt <- diag(-1, nrow = nL-1)
  mt <- -1*((row(mt) == col(mt) + (1)) + 0)
  diag(mt) <- rep(1, nL-1)
  et <- rep(linkf$dfun(eta), each=length(eta)) * mt
  et
}

pGumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE)
{
  q <- (q - loc)/scale
  p <- exp(-exp(q))
  if (lower.tail)
    1 - p
  else p
}

dGumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
  x <- -(x - loc)/scale
  d <- log(1/scale) - x - exp(-x)
  if (!log)
    exp(d)
  else d
}

qGumbel <- function (p, location = 0, scale = 1, lower.tail = TRUE, max = FALSE)
{
  if (!lower.tail)
    p <- 1 - p
  if (max)
    location - scale * log(-log(p))
  else location + scale * log(-log(1 - p))
}

yMx <- function(y, obs, nL)
{
  ym <- matrix(0, nrow = obs, ncol = nL, dimnames = list(NULL, levels(y)))
  yi <- as.integer(y)
  ym[cbind(1:obs, yi)] <- 1
  ym
}

lnkfun <- function(link){
  structure(list(
    pfun = switch(link, logit = plogis, probit = pnorm,
                  loglog = pgumbel, cloglog = pGumbel, cauchit = pcauchy),
    dfun = switch(link, logit = dlogis, probit = dnorm,
                  loglog = dgumbel, cloglog = dGumbel, cauchit = dcauchy),
    qfun = switch(link, logit = qlogis, probit = qnorm,
                  loglog = qgumbel, cloglog = qGumbel, cauchit = qcauchy)
  ), name = "linkfun")
}

CVserp <- function(lambda, dt, nFold, linkf, link, m,
                   slope, globalvar, nv, reverse, vnull, control, subs, wt)
{
  Indexes <- unlist(subs, use.names = FALSE)
  tnData  <- dt[-Indexes, ]
  tsData  <- dt[ Indexes, ]
  y <- as.ordered(tnData[, 1])
  x <- model.matrix(~tnData[, -1])
  if (is.null(dim(x))) x <- cbind(x)
  nL <- nL1 <- length(levels(y))
  ylev1 <- levels(y)[-length(levels(y))]
  obs <- length(y)
  yMtx <- yMx(y, obs, nL)
  yFreq <- colSums(yMtx)/obs
  xlst <- formxL(x, nL, slope, globalvar, m, vnull)
  xMat <- do.call(rbind, xlst)
  inibeta <- startv(link, linkf, globalvar, x, yFreq, xMat, nL, nv, slope)
  npar <- length(inibeta)
  x <- tnData[, -1L, drop = FALSE]
  wt <- rep(1, obs)
  estx <- serp.fit(lambda, globalvar, x, y, inibeta, xlst,
                   xMat, yMtx, nL, obs, npar, linkf, link,
                   reverse, vnull, control, slope, wt)
  est <- estx$coef
  if (!reverse) est <- -1L*est
  y2 <- as.ordered(tsData[, 1L])
  x2 <- model.matrix(~tsData[, -1L])
  nL <- length(levels(y2))
  obs <- length(y2)
  xlst <- formxL(x2, nL, slope, globalvar, m, vnull)
  xMat <- do.call(rbind, xlst)
  ylev2 <- levels(y2)[-length(levels(y2))]
  if (length(ylev1) != length(ylev2))
  {
    sk <- function(dff) seq(dff, length(est), (nL1 - 1))
    dff <- matrix(setdiff(ylev1, ylev2))
    est <- est[ -c(apply(dff, 1, sk))]
  }
  pr <- prlg(est, xMat, obs, yMtx = NULL, penx = NULL, linkf,
             control = NULL, wt = NULL)$pr
  double.eps <- sqrt(.Machine$double.xmin)
  pr[pr <= double.eps] <- double.eps
  yMtx <- yMx(y2, obs, nL)
  rs <- rowSums(yMtx)
  bs <- sum(yMtx * (1 - pr)^2 + (rs - yMtx) * pr^2) / sum(rs)
  bs
}

cvfun <- function(lambda, x, y, nFold, linkf, link, m,
                  slope, globalvar, nv, reverse, vnull, control, wt){
  df <- cbind(y, x)
  set.seed(111)
  dt <- df[sample(nrow(df)), ]
  obs <- dim(dt)[1L]
  folds  <- cut(seq(1, obs), breaks=nFold, labels=FALSE)
  fb <- data.frame(1:obs)
  sp <- split(fb, folds)
  tryCatch({
    rs <- mean(sapply(X = sp, FUN = CVserp, lambda = lambda,
                      dt = dt, nFold = nFold, linkf = linkf, link = link,
                      m = m, slope = slope, globalvar = globalvar, nv = nv,
                      reverse = reverse, vnull = vnull, control = control, wt = wt));},
    error=function(e) {rs <<- NA})
  rs
}

dvfun <- function(lambda, globalvar, x, y, inibeta,xlst, xMat,
                  yMtx, nL, obs, npar, linkf, link, reverse,
                  vnull,control, slope, wt, tuning){
  tryCatch({
    rr <- serp.fit(lambda, globalvar, x, y, inibeta,xlst, xMat,
                   yMtx, nL, obs, npar, linkf, link, reverse,
                   vnull,control, slope, wt);},
    error=function(e) {rr <<- NA})
  if(tuning == "finite")
    tryCatch({rd <- -2*as.numeric(rr$problik$llik2);}, error=function(e) {rd <<- NA})
  else
    tryCatch({rd <- rr$deviance;}, error=function(e) {rd <<- NA})
  rd

}
