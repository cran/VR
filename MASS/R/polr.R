# file MASS/polr.q
# copyright (C) 1994-2001 W. N. Venables and B. D. Ripley
#
polr <- function(formula, data = NULL, weights, start, ..., subset,
                 na.action = na.fail, contrasts = NULL, Hess=FALSE)
{
    logit <- function(p) log(p/(1-p))

    fmin <- function(beta) {
        gamm <- c(-100, beta[pc+1:q], 100)
        eta <- offset
        if(pc > 0) eta <- eta + drop(x %*% beta[1:pc])
        pr <- plogis(gamm[y+1] - eta) - plogis(gamm[y] - eta)
        if(all(pr > 0)) -sum(wt * log(pr)) else Inf
    }

    gmin <- function(beta) {
        gamm <- c(-100, beta[pc+1:q], 100)
        eta <- offset
        if(pc > 0) eta <- eta + drop(x %*% beta[1:pc])
        pr <- plogis(gamm[y+1] - eta) - plogis(gamm[y] - eta)
        p1 <- dlogis(gamm[y+1] - eta)
        p2 <- dlogis(gamm[y] - eta)
        g1 <- if(pc > 0) t(x) %*% (wt*(p1-p2)/pr) else numeric(0)
        xx <- .polrY1*p1 - .polrY2*p2
        d <- pmin(diff(beta[pc+1:q]), 0)
        g2 <- - t(xx) %*% (wt/pr)
        if(all(pr) > 0) c(g1, g2) else rep(NA, pc+q)
    }

    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    n <- nrow(x)
    pc <- ncol(x)
    if(xint > 0) {
        x <- x[, -xint, drop=FALSE]
        pc <- pc - 1
    } else warning("an intercept is needed and assumed")
    wt <- model.extract(m, weights)
    if(!length(wt)) wt <- rep(1, n)
    offset <- model.extract(m, offset)
    if(length(offset) <= 1) offset <- rep(0, n)
    y <- model.extract(m, response)
    if(!is.factor(y)) stop("response must be a factor")
    lev <- levels(y)
    if(length(lev) <= 2) stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    assign(".polrY1", col(Y) == y)
    assign(".polrY2", col(Y) == y-1)
    if(missing(start)) {
        # try logistic regression on `middle' cut
        q1 <- length(lev) %/% 2
        y1 <- (y > q1)
        X <- cbind(Intercept = rep(1, n), x)
        fit <- glm.fit(X, y1, wt, family = binomial(), offset = offset)
        coefs <- fit$coefficients
        spacing <- logit((1:q)/(q+1))
        start <- c(coefs[-1], -coefs[1] + spacing - spacing[q1])
    }
    res <- optim(start, fmin, gmin, method="BFGS", hessian = Hess)
    beta <- res$par[seq(len=pc)]
    zeta <- res$par[pc + 1:q]
    deviance <- 2 * res$value
    niter <- c(f.evals=res$counts[1], g.evals=res$counts[2])
    names(zeta) <- paste(lev[-length(lev)], lev[-1], sep="|")
    if(pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    } else {
        eta <- rep(0, n)
    }
    cumpr <- matrix(plogis(matrix(zeta, n, q, byrow=TRUE) - eta), , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
                fitted.values = fitted, lev = lev, terms = Terms,
                df.residual = sum(wt) - pc - q, edf = pc + q, n = sum(wt),
                call = match.call(), convergence = res$code, niter=niter)
    if(Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    attr(fit, "na.message") <- attr(m, "na.message")
    if(!is.null(attr(m, "na.action"))) fit$na.action <- attr(m, "na.action")
    fit$contrasts <- attr(x, "contrasts")
    fit$xlevels <- xlev
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
    cat("\nResidual Deviance:", format(round(x$deviance, 2)), "\n")
    cat("AIC:", format(round(x$deviance + 2*x$edf, 2)), "\n")
    invisible(x)
}

vcov.polr <- function(object)
{
    if(is.null(object$Hessian)) {
        cat("\nRe-fitting to get Hessian\n\n")
        object <- update(object, Hess=TRUE, trace=FALSE,
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
        print(x$coef[seq(len=pc), ], quote = FALSE, ...)
    } else {
        cat("\nNo coefficients\n")
    }
    cat("\nIntercepts:\n")
    print(coef[(pc+1):nrow(coef), ], quote = FALSE, ...)
    cat("\nResidual Deviance:", format(round(x$deviance, 2)), "\n")
    cat("AIC:", format(round(x$deviance + 2*x$edf, 2)), "\n")
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
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(X), nomatch=0)
        if(xint > 0) X <- X[, -xint, drop=FALSE]
        n <- nrow(X)
        q <- length(object$zeta)
        eta <- drop(X %*% object$coef)
        cumpr <- matrix(plogis(matrix(object$zeta, n, q, byrow=TRUE) - eta), , q)
        Y <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
        dimnames(Y) <- list(rownames(X), object$lev)
    }
    if(missing(newdata) && !is.null(object$na.action))
        Y <- napredict(object$na.action, Y)
    switch(type, class={
        Y <- factor(max.col(Y), levels=seq(along=object$lev),
                    labels=object$lev)
    }, probs={})
    drop(Y)
}

extractAIC.polr <- function(fit, scale = 0, k = 2, ...)
{
    edf <- fit$edf
    c(edf, deviance(fit) + k * edf)
}

model.frame.polr <- function(formula, data, na.action, ...)
{
    m <- formula$call
    m$start <- m$Hess <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    data <- eval(m, parent.frame())
    if(!is.null(mw <- m$weights)) {
        nm <- names(data)
        nm[match("(weights)", nm)] <- as.character(mw)
        names(data) <- nm
    }
    data
}
