# file MASS/vcov.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
vcov <- function(object, ...) UseMethod("vcov")
vcov.nls <- function(object)
{
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}
vcov.glm <- function(object)
{
    so <- summary(object, corr=FALSE)
    so$dispersion * so$cov.unscaled
}
vcov.lm <- function(object)
{
    so <- summary(object, corr=FALSE)
    so$sigma^2 * so$cov.unscaled
}

#deviance.nls <- function(object) sum(object$residuals^2)
deviance.nls <- function(object) object$m$deviance()

vcov.nlregb <- function(object, method=c("Fisher", "observed", "Huber"),
   scale=object$scale, eps=1e-3, tol=1)
{
    vcovscale <- function(devfn, par, scale, tol, ...)
    {
        # find a suitable initial scaling
        p <- length(par)
        ind <- 1:p
        v0 <- sum(devfn(par, ...)^2)
        scale <- sqrt(tol)/scale
        for (i in ind) {
            eps <- scale[i]
            inc <- ind == i
            v <- sum(devfn(par + eps*inc, ...)^2)
            if(v - v0 > tol) {
                repeat {
                    eps <- eps/3
                    v <- sum(devfn(par + eps*inc, ...)^2)
                    if(v - v0 < tol) break
                }
            } else {
                repeat {
                    eps <- eps*3
                    if(eps > 1000*scale[i])
                        stop(paste("scaling on var",i,"is too small"))
                    v <- sum(devfn(par + eps*inc, ...)^2)
                    if(v - v0 > tol) {eps <- eps/3; break}
                }
            }
            scale[i] <- eps
        }
        scale
    }

    method <- match.arg(method)
    par <- object$param
    n <- length(object$resid)
    p <- length(object$param)
    if(length(scale) == 1) scale <- rep(scale, p)
    s2 <- sum(object$resid^2)/(n-p)
    if(!is.null(gr <- object$jacobian)) {
        # gradient supplied
        K <- t(gr) %*% gr
        if(method == "Fisher")
            H2 <- 0
        else {
            r <- t(object$resid)
            grfn <- eval(substitute(object$call$jacobian))
            grfn <- eval(as.expression(grfn), sys.frame(sys.parent()))
            if(!length(grfn)) stop("Jacobian fn not found")
            sc <- eps * pmin(1, abs(par)) * sign(par)
            argnames <- names(grfn)
            argnames <- argnames[2:(length(argnames) - 1)]
            addargs <- object$aux[argnames]
            g <- object$jacobian
            H2 <- matrix(0, p, p)
            for(i in 1:p) {
                p1 <- par
                p1[i] <- p1[i] + sc[i]
                g1 <- do.call("grfn", c(list(p1), addargs))
                H2[i, ] <- r %*% (g1 - g)/sc[i]
            }
            H2 <- 0.5*(H2 + t(H2))
        }
        Jinv <- solve(K + H2)
        if(method != "Huber") V <- s2 * Jinv
        else V <- s2 * Jinv %*% K %*% Jinv
    } else {
        # no gradient supplied
        if(method != "Fisher") warning("Only Fisher information is available without gradient information")
        fn <- eval(substitute(object$call$residuals))
        fn <- eval(as.expression(fn), sys.frame(sys.parent()))
        if(!length(fn)) stop("objective fn not found")
        argnames <- names(fn)
        argnames <- argnames[-c(1, length(argnames))]
        addargs <- object$aux[argnames]
        scale <- do.call("vcovscale", c(list(fn, par, scale=scale, tol=tol*s2), addargs))
        ind <- 1:p
        H <- matrix(, n, p)
        for (j in 1:p)
            H[,j] <- do.call("fn", c(list(par + scale[j]*(ind==j)), addargs)) -
                do.call("fn", c(list(par - scale[j]*(ind==j)), addargs))
        J <- 0.25*crossprod(H)/outer(scale, scale)
        V <- s2 * solve(J)
    }
    v <- 2*sqrt(diag(V))
    upper <- eval(substitute(object$call$upper))
    if(is.null(upper)) upper <- Inf else upper <- eval(upper, sys.frame(sys.parent()))
    lower <- eval(substitute(object$call$lower))
    if(is.null(lower)) lower <- -Inf else lower <- eval(lower, sys.frame(sys.parent()))
    if(any(par - v <= lower) || any(par + v >= upper))
        warning("estimate is near the boundary: the estimated variance matrix may not be valid")
    V
}

vcov.nlminb <- function(object, tol=1, scale=object$scale, eps=1e-3, eps0=1)
{
    vcovscale <- function(devfn, par, scale, tol, ...)
    {
        # find a suitable initial scaling
        p <- length(par)
        ind <- 1:p
        v0 <- devfn(par, ...)
        scale <- sqrt(tol)/scale
        for (i in ind) {
            eps <- scale[i]
            inc <- ind == i
            v <- devfn(par + eps*inc, ...)
            if(v - v0 > tol) {
                repeat {
                    eps <- eps/3
                    v <- devfn(par + eps*inc, ...)
                    if(v - v0 < tol) break
                }
            } else {
                repeat {
                    eps <- eps*3
                    if(eps > 1000*scale[i])
                        stop(paste("scaling on var",i,"is too small"))
                    v <- devfn(par + eps*inc, ...)
                    if(v - v0 > tol) {eps <- eps/3; break}
                }
            }
            scale[i] <- eps
        }
        scale
    }

    vcov0 <- function(devfn, par, scale = rep(1, p), ...)
    {
        p <- length(par)
        ind <- 1:p
        v0 <- devfn(par, ...)
        h <- numeric(p)
        H <- matrix(, p, p)
        for (i in 1:p)
            h[i] <- devfn(par + scale[i]*(ind==i), ...) - v0
        if(p > 1) for (i in 2:p) {
            inc <- scale[i] * (ind == i)
            for(j in 1:(i-1))
                H[i, j] <- H[j, i] <-
                    devfn(par + inc + scale[j]*(ind==j), ...) - v0
        }
        H <- H - outer(h, h, "+")
        diag(H) <- 2*h
        H/outer(scale, scale)
    }

    par <- object$param
    p <- length(par)
    fn <- eval(substitute(object$call$objective))
    fn <- eval(as.expression(fn), sys.frame(sys.parent()))
    if(!length(fn)) stop("objective fn not found")
    argnames <- names(fn)
    argnames <- argnames[-c(1, length(argnames))]
    addargs <- object$aux[argnames]
    if(length(scale) == 1) scale <- rep(scale, p)
    hessian <- object$call$hessian
    # object$hessian appears to lie, so we redo the calculation
    if(is.logical(hessian) && hessian) {
        # hessian in gradient
        grfn <- eval(substitute(object$call$gradient))
        grfn <- eval(as.expression(grfn), sys.frame(sys.parent()))
        if(!length(grfn)) stop("gradient fn not found")
        hh <- do.call("grfn", c(list(par), addargs))$hessian
        H <- matrix(0, p, p)
        H[row(H) <= col(H)] <- hh
        H[row(H) > col(H)] <- t(H)[row(H) > col(H)]
    } else if(!is.logical(hessian) && !is.null(hessian)) {
        # separate Hessian function
        hfn <- eval(substitute(object$call$hessian))
        hfn <- eval(as.expression(hfn), sys.frame(sys.parent()))
        if(!length(hfn)) stop("hessian fn not found")
        hh <- do.call("hfn", c(list(par), addargs))
        H <- matrix(0, p, p)
        H[row(H) <= col(H)] <- hh
        H[row(H) > col(H)] <- t(H)[row(H) > col(H)]
    } else {
        scale <- do.call("vcovscale", c(list(fn, par, scale=scale, tol=tol), addargs))
        if(!is.null(object$call$gradient)) {
            # gradient but no Hessian
            grfn <- eval(substitute(object$call$gradient))
            grfn <- eval(as.expression(grfn), sys.frame(sys.parent()))
            if(!length(grfn)) stop("gradient fn not found")
            sc <- eps * scale
            g <- do.call("grfn", c(list(par), addargs))
            H <- matrix(0, p, p)
            for(i in 1:p) {
                g1 <- do.call("grfn", c(list(par - sc[i]*(1:p==i)), addargs))
                g2 <- do.call("grfn", c(list(par + sc[i]*(1:p==i)), addargs))
                H[i, ] <- 0.5*(g2 - g1)/sc[i]
            }
        } else {
            # No derivatives available
            H <- do.call("vcov0", c(list(fn, par, scale=eps0*scale), addargs))
        }
    }
    V <- solve(0.5*(H + t(H)))
    if(any(diag(V) < 0)) stop("Estimated variances are negative")
    v <- 2*sqrt(diag(V))
    upper <- eval(substitute(object$call$upper))
    if(is.null(upper)) upper <- Inf else upper <- eval(upper, sys.frame(sys.parent()))
    lower <- eval(substitute(object$call$lower))
    if(is.null(lower)) lower <- -Inf else lower <- eval(lower, sys.frame(sys.parent()))
    if(any(par - v <= lower) || any(par + v >= upper))
        warning("estimate is near the boundary: the estimated variance matrix may not be valid")
    V
}
