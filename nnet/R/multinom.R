# file nnet/multinom.q copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#

multinom <-
    function(formula, data, weights, subset, na.action,
             contrasts = NULL, Hess = FALSE, summ = 0, censored = FALSE,
             model = FALSE, ...)
{
  class.ind <- function(cl)
  {
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (as.vector(unclass(cl)) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  summ2 <- function(X, Y)
  {
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Y)
    Z <- t(cbind(X, Y))
    storage.mode(Z) <- "double"
    z <- .C("VR_summ2",
	    as.integer(n),
	    as.integer(p),
	    as.integer(q),
	    Z = Z,
	    na = integer(1), PACKAGE = "nnet")
    Za <- t(z$Z[, 1:z$na, drop = FALSE])
    list(X = Za[, 1:p, drop = FALSE], Y = Za[, p + 1:q])
  }

  call <- match.call()
  m <- match.call(expand = FALSE)
  m$summ <- m$Hess <- m$contrasts <- m$censored <- m$model <- m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  X <- model.matrix(Terms, m, contrasts)
  cons <- attr(X, "contrasts")
  Xr <- qr(X)$rank
  Y <- model.response(m)
  if(!is.matrix(Y)) Y <- as.factor(Y)
  w <- model.weights(m)
  if(length(w) == 0)
    if(is.matrix(Y)) w <- rep(1, dim(Y)[1])
    else w <- rep(1, length(Y))
  lev <- levels(Y)
  if(is.factor(Y)) {
    counts <- table(Y)
    if(any(counts == 0)) {
      warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),
		    "are empty"))
      Y <- factor(Y, levels=lev[counts > 0])
      lev <- lev[counts > 0]
    }
    if(length(lev) == 2) Y <- as.vector(unclass(Y)) - 1
    else Y <- class.ind(Y)
  }
  if(summ==1) {
    Z <- cbind(X, Y)
    assign("z1", cumprod(apply(Z, 2, max)+1))
    Z1 <- apply(Z, 1, function(x) sum(z1*x))
    oZ <- order(Z1)
    Z2 <- !duplicated(Z1[oZ])
    oX <- (seq(along=Z1)[oZ])[Z2]
    X <- X[oX, , drop=FALSE]
    Y <- if(is.matrix(Y)) Y[oX, , drop=FALSE] else Y[oX]
    w <- diff(c(0,cumsum(w))[c(Z2,TRUE)])
    print(dim(X))
  }
  if(summ==2) {
    Z <- summ2(cbind(X, Y), w)
    X <- Z$X[, 1:ncol(X)]
    Y <- Z$X[, ncol(X) + 1:ncol(Y), drop = FALSE]
    w <- Z$Y
    print(dim(X))
  }
  if(summ==3) {
    Z <- summ2(X, Y*w)
    X <- Z$X
    Y <- Z$Y[, 1:ncol(Y), drop = FALSE]
    w <- rep(1, nrow(X))
    print(dim(X))
  }
  offset <- model.offset(m)
  r <- ncol(X)
  if(is.matrix(Y)) {
# 3 or more response levels or direct matrix spec.
    p <- ncol(Y)
    sY <- Y %*% rep(1, p)
    if(any(sY==0)) stop("some case has no observations")
    if(!censored) {
      Y <- Y / matrix(sY, nrow(Y), p)
      w <- w*sY
    }
    if(length(offset) > 1) {
      if(ncol(offset) !=  p) stop("ncol(offset) is wrong")
      mask <- c(rep(0, r+1+p), rep(c(0, rep(1, r), rep(0, p)), p-1) )
      X <- cbind(X, offset)
      Wts <- as.vector(rbind(matrix(0, r+1, p), diag(p)))
      fit <- nnet.default(X, Y, w, Wts=Wts, mask=mask, size=0, skip=TRUE,
                          softmax=TRUE, censored=censored, rang=0, ...)
    } else {
      mask <- c(rep(0, r+1), rep(c(0, rep(1, r)), p-1) )
      fit <- nnet.default(X, Y, w, mask=mask, size=0, skip=TRUE, softmax=TRUE,
                          censored=censored, rang=0, ...)
    }
  } else {
# 2 response levels
    if(length(offset) <= 1) {
      mask <- c(0, rep(1, r))
      fit <- nnet.default(X, Y, w, mask=mask, size=0, skip=TRUE, entropy=TRUE,
                          rang=0, ...)
    } else {
      mask <- c(0, rep(1, r), 0)
      Wts <- c(rep(0, r+1), 1)
      X <- cbind(X, offset)
      fit <- nnet.default(X, Y, w, Wts=Wts, mask=mask, size=0, skip=TRUE,
                  entropy=TRUE, rang=0, ...)
    }
  }
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$terms <- Terms
  fit$call <- call
  fit$weights <- w
  fit$lev <- lev
  fit$deviance <- 2 * fit$value
  fit$rank <- Xr
  edf <- ifelse(length(lev) == 2, 1, length(lev)-1)*Xr
  if(is.matrix(Y)) {
    edf <- (ncol(Y)-1)*Xr
    if(length(dn <- colnames(Y)) > 0) fit$lab <- dn
    else fit$lab <- 1:ncol(Y)
  }
  fit$coefnames <- colnames(X)
  fit$vcoefnames <- fit$coefnames[1:r] # remove offset cols
  fit$na.action <- attr(m, "na.action")
  fit$contrasts <- cons
  fit$xlevels <- .getXlevels(Terms, m)
  fit$edf <- edf
  fit$AIC <- fit$deviance + 2 * edf
  if(model) fit$model <- m
  class(fit) <- c("multinom", "nnet")
  if(Hess) fit$Hessian <- multinomHess(fit, X)
  fit
}

predict.multinom <- function(object, newdata, type=c("class","probs"), ...)
{
  if(!inherits(object, "multinom")) stop("Not a multinom fit")
  type <- match.arg(type)
  if(missing(newdata)) Y <- fitted(object)
  else {
    newdata <- as.data.frame(newdata)
    rn <- row.names(newdata)
    Terms <- delete.response(object$terms)
    m <- model.frame(Terms, newdata, na.action = na.omit,
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")) &&
        exists(".checkMFClasses", envir=NULL)) .checkMFClasses(cl, m)
    keep <- match(row.names(m), rn)
    X <- model.matrix(Terms, m, contrasts = object$contrasts)
    Y1 <- predict.nnet(object, X)
    Y <- matrix(NA, nrow(newdata), ncol(Y1),
                dimnames = list(rn, colnames(Y1)))
    Y[keep, ] <- Y1
  }
  switch(type, class={
    if(length(object$lev) > 2)
      Y <- factor(max.col(Y), levels=seq(along=object$lev), labels=object$lev)
    if(length(object$lev) == 2)
      Y <- factor(1 + (Y > 0.5), levels=1:2, labels=object$lev)
    if(length(object$lev) == 0)
      Y <- factor(max.col(Y), levels=seq(along=object$lab), labels=object$lab)
  }, probs={})
  drop(Y)
}

print.multinom <- function(x, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nCoefficients:\n")
  print(coef(x), ...)
  cat("\nResidual Deviance:", format(x$deviance), "\n")
  cat("AIC:", format(x$AIC), "\n")
  invisible(x)
}

coef.multinom <- function(object, ...)
{
  r <- length(object$vcoefnames)
  if(length(object$lev) == 2) {
    coef <- object$wts[1+(1:r)]
    names(coef) <- object$vcoefnames
  } else {
    coef <- matrix(object$wts, nrow = object$n[3], byrow=TRUE)[, 1+(1:r), drop=FALSE]
    if(length(object$lev)) dimnames(coef) <- list(object$lev, object$vcoefnames)
    if(length(object$lab)) dimnames(coef) <- list(object$lab, object$vcoefnames)
    coef <- coef[-1, , drop=FALSE]
  }
  coef
}

drop1.multinom <- function(object, scope, sorted = FALSE, trace = FALSE, ...)
{
  if(!inherits(object, "multinom")) stop("Not a multinom fit")
  if(missing(scope)) scope <- drop.scope(object)
    else {
      if(!is.character(scope))
	scope <- attr(terms(update.formula(object, scope)), "term.labels")
      if(!all(match(scope, attr(object$terms, "term.labels"),
                    nomatch = FALSE)))
	stop("scope is not a subset of term labels")
    }
  ns <- length(scope)
  ans <- matrix(nrow = ns+1, ncol = 2,
                dimnames = list(c("<none>", scope), c("Df", "AIC")))
  ans[1, ] <- c(object$edf, object$AIC)
  i <- 2
  for(tt in scope) {
    cat("trying -", tt,"\n")
    nobject <- update(object, paste("~ . -", tt), trace = trace)
    if(nobject$edf == object$edf) nobject$AIC <- NA
    ans[i, ] <- c(nobject$edf, nobject$AIC)
    i <- i+1
  }
  if(sorted) ans <- ans[order(ans[, 2]), ]
  as.data.frame(ans)
}

add1.multinom <- function(object, scope, sorted = FALSE, trace = FALSE, ...)
{
  if(!inherits(object, "multinom")) stop("Not a multinom fit")
  if(!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope,
					   evaluate = FALSE))
  if(!length(scope))
    stop("no terms in scope for adding to object")
  ns <- length(scope)
  ans <- matrix(nrow = ns+1, ncol = 2,
                dimnames = list(c("<none>",paste("+",scope,sep="")),
                  c("Df", "AIC")))
  ans[1, ] <- c(object$edf, object$AIC)
  i <- 2
  for(tt in scope) {
    cat("trying +", tt,"\n")
    nobject <- update(object, paste("~ . +", tt), trace=trace)
    if(nobject$edf == object$edf) nobject$AIC <- NA
    ans[i, ] <- c(nobject$edf, nobject$AIC)
    i <- i+1
  }
  if(sorted) ans <- ans[order(ans[, 2]), ]
  as.data.frame(ans)
}

extractAIC.multinom <- function(fit, scale, k = 2, ...)
  c(fit$edf, fit$AIC + (k-2)*fit$edf)

vcov.multinom <- function(object, ...)
{
  ginv <- function(X, tol = sqrt(.Machine$double.eps))
    {
    #
    # simplified version of ginv in MASS
    #
      Xsvd <- svd(X)
      Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
      if(!any(Positive)) array(0, dim(X)[2:1])
      else Xsvd$v[, Positive] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive]))
    }

  if(is.null(Hess <- object$Hessian)) Hess <- multinomHess(object)
  structure(ginv(Hess), dimnames = dimnames(Hess))
}

summary.multinom <-
function(object, correlation = TRUE, digits = options()$digits,
         Wald.ratios = FALSE, ...)
{
  vc <- vcov(object)
  r <- length(object$vcoefnames)
  se <- sqrt(diag(vc))
  if(length(object$lev) == 2) {
    coef <- object$wts[1 + (1:r)]
    stderr <- se
    names(coef) <- names(stderr) <- object$vcoefnames
  } else {
    coef <- matrix(object$wts, nrow = object$n[3],
		   byrow = TRUE)[-1, 1 + (1:r), drop = FALSE]
    stderr <- matrix(se, nrow = object$n[3] - 1, byrow = TRUE)
    if(length(l <- object$lab) || length(l <- object$lev))
      dimnames(coef) <- dimnames(stderr) <- list(l[-1], object$vcoefnames)
  }
  object$is.binomial <- (length(object$lev) == 2)
  object$digits <- digits
  object$coefficients <- coef
  object$standard.errors <- stderr
  if(Wald.ratios) object$Wald.ratios <- coef/stderr
  if(correlation) object$correlation <- vc/outer(se, se)
  class(object) <- "summary.multinom"
  object
}

print.summary.multinom <- function(x, digits = x$digits, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  cat("\nCoefficients:\n")
  if(x$is.binomial) {
    print(cbind(Values = x$coefficients,
		"Std. Err." = x$standard.errors,
		"Value/SE" = x$Wald.ratios),
	  digits = digits)
  } else {
    print(x$coefficients, digits = digits)
    cat("\nStd. Errors:\n")
    print(x$standard.errors, digits = digits)
    if(!is.null(x$Wald.ratios)) {
      cat("\nValue/SE (Wald statistics):\n")
      print(x$coefficients/x$standard.errors, digits = digits)
    }
  }
  cat("\nResidual Deviance:", format(x$deviance), "\n")
  cat("AIC:", format(x$AIC), "\n")
  if(!is.null(correl <- x$correlation)) {
    p <- dim(correl)[2]
    if(p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      ll <- lower.tri(correl)
      correl[ll] <- format(round(correl[ll], digits))
      correl[!ll] <- ""
      print(correl[-1, -p], quote = FALSE, ...)
    }
  }
  invisible(x)
}

anova.multinom <- function(object, ..., test = c("Chisq", "none"))
{
  test <- match.arg(test)
  dots <- list(...)
  if(length(dots) == 0)
    stop("anova is not implemented for a single multinom object")
  mlist <- list(object, ...)
  nt <- length(mlist)
  dflis <- sapply(mlist, function(x) x$edf)
  s <- order(dflis)
  dflis <- nrow(residuals(object)) * (ncol(residuals(object))-1) - dflis
  mlist <- mlist[s]
  if(any(!sapply(mlist, inherits, "multinom")))
    stop("not all objects are of class 'multinom'")
  rsp <- unique(sapply(mlist, function(x) paste(formula(x)[2])))
  mds <- sapply(mlist, function(x) paste(formula(x)[3]))
  dfs <- dflis[s]
  lls <- sapply(mlist, function(x) deviance(x))
  tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
  out <- data.frame(Model = mds, Resid.df = dfs,
                    Deviance = lls, Test = tss, Df = df, LRtest = x2,
                    Prob = pr)
  names(out) <- c("Model", "Resid. df", "Resid. Dev", "Test",
                  "   Df", "LR stat.", "Pr(Chi)")
  if(test=="none") out <- out[, 1:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of Multinomial Models\n",
      paste("Response:", rsp))
  out
}


model.frame.multinom <- function(formula, ...)
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    if(any(nargs > 0) || is.null(formula$model)) {
        oc <- formula$call
        oc[[1]] <- as.name("model.frame")
        m <- match(names(oc)[-1], c("formula", "data", "na.action", "subset"))
        oc <- oc[c(TRUE, !is.na(m))]
        oc[names(nargs)] <- nargs
        if (is.null(env <- environment(formula$terms))) env <- parent.frame()
        eval(oc, env)
    } else formula$model
}
