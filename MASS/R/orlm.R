# file MASS/orlm.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
orlm <-
    function(formula, data, weights, subset, na.action, method = "qr",
             model = FALSE, x = FALSE, y = FALSE, ...)
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$method <- m$model <- m$x <- m$y <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    if(method == "model.frame") return(m)
    Terms <- attr(m, "terms")
    weights <- model.extract(m, weights)
    Y <- model.extract(m, response)
    X <- model.matrix(Terms, m)
    if(length(weights)==0) weights <- rep(1, nrow(X))
    fit <- hsreg(X, Y, wx = weights, ...)
    fit$terms <- Terms
    fit$call <- call
    if(model) fit$model <- m
    fit$x <- X
    fit$y <- Y
    attr(fit, "na.message") <- attr(m, "na.message")
    if(!is.null(attr(m, "na.action"))) fit$na.action <- attr(m, "na.action")
    fit
}

hsreg <-
    function(x, y, w = rep(1, nrow(x)), k=1.345, wx, maxit = 20, sw=1000,
             acc = .Machine$double.eps^0.25, test.vec = "resid", ...)
{
    irls.delta <- function(old, new)
        sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
    irls.rrxwr <- function(x, w, r)
    {
        w <- sqrt(w)
        max(abs((matrix(r * w, 1, length(r)) %*% x)/
                sqrt(matrix(w, 1, length(r)) %*% (x^2))))/sqrt(sum(w * r^2))
    }
    if(!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec)))
        stop("invalid testvec")
    if(!missing(wx)) {
        if(length(wx) != nrow(x))
            stop("Length of wx must equal number of observations")
        if(any(wx < 0))
            stop("Negative wx value")
        w <- w * wx
    }
    temp <- lm.wfit(x, y, w, method="qr", ...)
    coef <- temp$coef
    resid <- temp$resid
    th <- 2*pnorm(k)-1
    gamma <- th + k^2*(1-th) -2*k*dnorm(k)

    done <- FALSE
    conv <- NULL
    n1 <- nrow(x) - ncol(x)
    scale <- median(abs(resid))/0.6745
    for(iiter in 1:maxit) {
        if(!is.null(test.vec)) testpv <- get(test.vec)
        ks <- k*scale
        if(iiter < sw) scale <- median(abs(resid))/0.6745
        else scale <- sqrt(sum(pmin(resid^2,ks^2))/(n1*gamma))
        if(scale == 0) {
            done <- TRUE
            break
        }
        w <- as.vector(wt.huber(resid/scale, k))
        if(!missing(wx)) w <- w * wx    # adjust for wx weights
        temp <- lm.wfit(x, y, w, method="qr")
        coef <- temp$coef
        resid <- temp$residuals
        if(!is.null(test.vec))
            convi <- irls.delta(testpv, get(test.vec))
        else convi <- irls.rrxwr(x, w, resid)
        conv <- c(conv, convi)
        done <- convi <= acc
        if(!done) next
        if(!exists("method.done") || method.done) break
    }
    if(!done)
        warning(paste("hsreg failed to converge in", maxit, "steps"))
    if(!missing(wx)) {
        tmp <- (wx != 0)
        w[tmp] <- w[tmp]/wx[tmp]
    }
    names(scale) <- NULL                # since median assigns a name
    fit <- list(coefficients = coef, residuals = resid,
                fitted.values = temp$fitted.values, rank = temp$rank,
                assign =temp$assign,  w = w, k = k, s = scale,
                conv = conv, converged = done)
    class(fit) <- c("orlm", "lm")
    fit
}

print.orlm <-
function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    if(x$converged)
        cat("Converged in",length(x$conv), "iterations\n")
    else cat("Ran",length(x$conv), "iterations without convergence\n")
    coef <- x$coef
    cat("\nCoefficients:\n")
    print(coef, ...)
    nobs <- length(x$resid)
    rdf <- nobs - length(coef)
    cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
    cat("Scale estimate:", format(signif(x$s,3)), "\n")
    invisible(x)
}

summary.orlm <- function(object, ...)
{
    k <- object$k
    s <- object$s
    ks <- k*s
    coef <- object$coef
    cnames <- names(coef)
    ptotal <- length(coef)
    resid <- object$resid
    n <- length(resid)
    if(any(na <- is.na(coef))) coef <- coef[!na]
    p <- length(coef)
    rdf <- n - p
    rinv <- diag(p)
    dimnames(rinv) <- list(cnames, cnames)
    S <- sum(pmin(abs(resid), ks)^2)/rdf
    m <- sum(abs(resid) < ks)/n
    kappa <- 1 + p*(1-m)/(n*m)
    stddev <- sqrt(S)*(kappa/m)
    R <- qr(object$x)$qr
    R <- R[1:p, 1:p, drop = FALSE]
    # for(i in 2:p)for(j in 1:(i-1))R[i,j] <- 0
    R[lower.tri(R)] <- 0
    rinv <- solve(R, rinv)
    rowlen <- (rinv^2 %*% rep(1, p))^0.5
    names(rowlen) <- cnames
    correl <- rinv * array(1/rowlen, c(p, p))
    correl <- correl %*% t(correl)
    coef <- array(coef, c(p, 3))
    dimnames(coef) <- list(cnames, c("Value", "Std. Error", "t value"))
    coef[, 2] <- rowlen %o% stddev
    coef[, 3] <- coef[, 1]/coef[, 2]
    object <- object["call"]
    object$residuals <- resid
    object$coefficients <- coef
    object$sigma <- s
    object$stddev <- stddev
    object$df <- c(p, rdf, ptotal)
    object$r.squared <- NA
    object$cov.unscaled <- rinv %*% t(rinv)
    object$correlation <- correl
    object$terms <- NA
    class(object) <- "summary.lm"
    object
}
