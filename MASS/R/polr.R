# file MASS/polr.q
# copyright (C) 1994-2004 W. N. Venables and B. D. Ripley
#
polr <- function(formula, data, weights, start, ..., subset,
                 na.action, contrasts = NULL, Hess = FALSE,
                 model = TRUE,
                 method = c("logistic", "probit", "cloglog", "cauchit"))
{
    logit <- function(p) log(p/(1 - p))

    fmin <- function(beta) {
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
        eta <- offset
        if (pc > 0)
            eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
        if (all(pr > 0))
            -sum(wt * log(pr))
        else Inf
    }
    jacobian <- function(theta) { ## dgamma by dtheta matrix
        k <- length(theta)
        etheta <- exp(theta)
        mat <- matrix(0 , k, k)
        mat[, 1] <- rep(1, k)
        for (i in 2:k) mat[i:k, i] <- etheta[i]
        mat
    }
    gmin <- function(beta) {
        theta <- beta[pc + 1:q]
        gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
        eta <- offset
        if(pc > 0) eta <- eta + drop(x %*% beta[1:pc])
        pr <- pfun(gamm[y+1] - eta) - pfun(gamm[y] - eta)
        p1 <- dfun(gamm[y+1] - eta)
        p2 <- dfun(gamm[y] - eta)
        g1 <- if(pc > 0) t(x) %*% (wt*(p1 - p2)/pr) else numeric(0)
        xx <- .polrY1*p1 - .polrY2*p2
        g2 <- - t(xx) %*% (wt/pr)
        g2 <- t(g2) %*% jacobian(theta)
        if(all(pr) > 0) c(g1, g2) else rep(NA, pc+q)
    }

    m <- match.call(expand.dots = FALSE)
    method <- match.arg(method)
    pfun <- switch(method, logistic = plogis, probit = pnorm,
                   cloglog = pgumbel, cauchit = pcauchy)
    dfun <- switch(method, logistic = dlogis, probit = dnorm,
                   cloglog = dgumbel, cauchit = dcauchy)
    if(is.matrix(eval.parent(m$data)))
        m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$method <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts") # will get dropped by subsetting
    if(xint > 0) {
        x <- x[, -xint, drop=FALSE]
        pc <- pc - 1
    } else warning("an intercept is needed and assumed")
    wt <- model.weights(m)
    if(!length(wt)) wt <- rep(1, n)
    offset <- model.offset(m)
    if(length(offset) <= 1) offset <- rep(0, n)
    y <- model.response(m)
    if(!is.factor(y)) stop("response must be a factor")
    lev <- levels(y)
    if(length(lev) <= 2) stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    if(missing(start)) {
        # try logistic/probit regression on 'middle' cut
        q1 <- length(lev) %/% 2
        y1 <- (y > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <-
            switch(method,
                   "logistic"= glm.fit(X, y1, wt, family = binomial(), offset = offset),
                   "probit" = glm.fit(X, y1, wt, family = binomial(probit), offset = offset),
                   ## this is deliberate, a better starting point
                   "cloglog" = glm.fit(X, y1, wt, family = binomial(probit), offset = offset),
                   "cauchit" = glm.fit(X, y1, wt, family = binomial(cauchit), offset = offset))
        if(!fit$converged)
            stop("attempt for find suitable starting values failed")
        coefs <- fit$coefficients
        if(any(is.na(coefs))) {
            warning("design appears to be rank-deficient, so dropping some coefs")
            keep <- names(coefs)[!is.na(coefs)]
            coefs <- coefs[keep]
            x <- x[, keep[-1], drop = FALSE]
            pc <- ncol(x)
        }
        spacing <- logit((1:q)/(q+1)) # just a guess
        if(method != "logit") spacing <- spacing/1.7
        gammas <- -coefs[1] + spacing - spacing[q1]
        thetas <- c(gammas[1], log(diff(gammas)))
        start <- c(coefs[-1], thetas)
    }
    res <- optim(start, fmin, gmin, method="BFGS", hessian = Hess, ...)
    beta <- res$par[seq(len=pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1],exp(theta[-1])))
    deviance <- 2 * res$value
    niter <- c(f.evals=res$counts[1], g.evals=res$counts[2])
    names(zeta) <- paste(lev[-length(lev)], lev[-1], sep="|")
    if(pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    } else {
        eta <- rep(0, n)
    }
    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow=TRUE) - eta), , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
                fitted.values = fitted, lev = lev, terms = Terms,
                df.residual = sum(wt) - pc - q, edf = pc + q, n = sum(wt),
                call = match.call(), method = method,
		convergence = res$convergence, niter = niter)
    if(Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if(model) fit$model <- m
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    class(fit) <- "polr"
    fit
}

print.polr <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    if(length(coef(x))) {
        cat("\nCoefficients:\n")
        print(coef(x), ...)
    } else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts:\n")
    print(x$zeta, ...)
    cat("\nResidual Deviance:", format(x$deviance, nsmall=2), "\n")
    cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2), "\n")
    if(x$convergence > 0)
        cat("Warning: did not converge as iteration limit reached\n")
    invisible(x)
}

vcov.polr <- function(object, ...)
{
    if(is.null(object$Hessian)) {
        cat("\nRe-fitting to get Hessian\n\n")
        object <- update(object, Hess=TRUE,
                         start=c(object$coef, object$zeta))
    }
    structure(ginv(object$Hessian), dimnames = dimnames(object$Hessian))

}

summary.polr <- function(object, digits = max(3, .Options$digits - 3),
                         correlation = FALSE, ...)
{
    cc <- c(coef(object), object$zeta)
    pc <- length(coef(object))
    q <- length(object$zeta)
    coef <- matrix(0, pc+q, 3, dimnames=list(names(cc),
                               c("Value", "Std. Error", "t value")))
    coef[, 1] <- cc
    vc <- vcov(object)
    z.ind <- (pc + 1):(pc + q)
    gamma <- object$zeta
    theta <- c(gamma[1], log(diff(gamma)))
    jacobian <- function(theta) { ## dgamma by dtheta matrix
        k <- length(theta)
        etheta <- exp(theta)
        mat <- matrix(0 , k, k)
        mat[, 1] <- rep(1, k)
        for (i in 2:k) mat[i:k, i] <- etheta[i]
        mat
    }
    J <- jacobian(theta)
    vc[z.ind, z.ind] <- J %*% vc[z.ind, z.ind] %*% t(J)
    coef[, 2] <- sd <- sqrt(diag(vc))
    coef[, 3] <- coef[, 1]/coef[, 2]
    object$coefficients <- coef
    object$pc <- pc
    object$digits <- digits
    if(correlation)
        object$correlation <- (vc/sd)/rep(sd, rep(pc+q, pc+q))
    class(object) <- "summary.polr"
    object
}

print.summary.polr <- function(x, digits = x$digits, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    coef <- format(round(x$coef, digits=digits))
    pc <- x$pc
    if(pc > 0) {
        cat("\nCoefficients:\n")
        print(x$coef[seq(len=pc), , drop=FALSE], quote = FALSE, ...)
    } else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts:\n")
    print(coef[(pc+1):nrow(coef), , drop=FALSE], quote = FALSE, ...)
    cat("\nResidual Deviance:", format(x$deviance, nsmall=2), "\n")
    cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2), "\n")
    if(!is.null(correl <- x$correlation)) {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        print(correl[-1, -ncol(correl)], quote = FALSE, ...)
    }
    invisible(x)
}

predict.polr <- function(object, newdata, type=c("class","probs"), ...)
{
    if(!inherits(object, "polr")) stop("Not a polr fit")
    type <- match.arg(type)
    if(missing(newdata)) Y <- object$fitted
    else {
        newdata <- as.data.frame(newdata)
        Terms <- delete.response(object$terms)
        m <- model.frame(Terms, newdata, na.action = function(x) x,
                         xlev = object$xlevels)
        if (!is.null(cl <- attr(Terms, "dataClasses")) &&
            exists(".checkMFClasses", envir=NULL)) .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(X), nomatch=0)
        if(xint > 0) X <- X[, -xint, drop=FALSE]
        n <- nrow(X)
        q <- length(object$zeta)
        eta <- drop(X %*% object$coef)
        pfun <- switch(object$method, logistic = plogis, probit = pnorm,
                       cloglog = pgumbel, cauchit = pcauchy)
        cumpr <- matrix(pfun(matrix(object$zeta, n, q, byrow=TRUE) - eta), , q)
        Y <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
        dimnames(Y) <- list(rownames(X), object$lev)
    }
    if(missing(newdata) && !is.null(object$na.action))
        Y <- napredict(object$na.action, Y)
    switch(type, class={
        Y <- factor(max.col(Y), levels=seq(along=object$lev),
                    labels=object$lev)
    }, probs = {})
    drop(Y)
}

extractAIC.polr <- function(fit, scale = 0, k = 2, ...)
{
    edf <- fit$edf
    c(edf, deviance(fit) + k * edf)
}

model.frame.polr <- function(formula, ...)
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
    if(length(nargs) || is.null(formula$model)) {
        m <- formula$call
        m$start <- m$Hess <- m$... <- NULL
        m[[1]] <- as.name("model.frame")
        m[names(nargs)] <- nargs
        if (is.null(env <- environment(formula$terms))) env <- parent.frame()
        data <- eval(m, env)
        if(!is.null(mw <- m$weights)) {
            nm <- names(data)
            nm[match("(weights)", nm)] <- as.character(mw)
            names(data) <- nm
        }
        data
    } else formula$model
}

pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
    q <- (q - loc)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 1 - p else p
}

dgumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
    d <- log(1/scale) - x - exp(-x)
    if (!log) exp(d) else d
}

anova.polr <- function (object, ..., test = c("Chisq", "none"))
{
    test <- match.arg(test)
    dots <- list(...)
    if (length(dots) == 0)
        stop("anova is not implemented for a single polr object")
    mlist <- list(object, ...)
    nt <- length(mlist)
    dflis <- sapply(mlist, function(x) x$edf)
    s <- order(dflis)
    mlist <- mlist[s]
    if (any(!sapply(mlist, inherits, "polr")))
        stop("not all objects are of class 'polr'")
    rsp <- unique(sapply(mlist, function(x) paste(formula(x)[2])))
    mds <- sapply(mlist, function(x) paste(formula(x)[3]))
    dfs <- dflis[s]
    lls <- sapply(mlist, function(x) deviance(x))
    tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
    df <- c(NA, diff(dfs))
    x2 <- c(NA, -diff(lls))
    pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
    out <- data.frame(Model = mds, Resid.df = dfs, Deviance = lls,
                      Test = tss, Df = df, LRtest = x2, Prob = pr)
    names(out) <- c("Model", "Resid. df", "Resid. Dev", "Test",
                    "   Df", "LR stat.", "Pr(Chi)")
    if (test == "none") out <- out[, 1:6]
    class(out) <- c("Anova", "data.frame")
    attr(out, "heading") <-
        c("Likelihood ratio tests of proportional odds models\n",
          paste("Response:", rsp))
    out
}
