# file MASS/rlm.q
# copyright (C) 1994-2004 W. N. Venables and B. D. Ripley
#
rlm <- function(x, ...) UseMethod("rlm")

rlm.formula <-
    function(formula, data, weights, ..., subset, na.action,
             method = c("M", "MM", "model.frame"),
             wt.method = c("inv.var", "case"),
             model = TRUE, x.ret = TRUE, y.ret = FALSE, contrasts = NULL)
{
    mf <- match.call(expand.dots = FALSE)
    mf$method <- mf$wt.method <- mf$model <- mf$x.ret <- mf$y.ret <- mf$contrasts <- mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    method <- match.arg(method)
    wt.method <- match.arg(wt.method)
    if(method == "model.frame") return(mf)
    mt <- attr(mf, "terms")
    y <- model.response(mf)
    x <- model.matrix(mt, mf, contrasts)
    xvars <- as.character(attr(mt, "variables"))[-1]
    if ((yvar <- attr(mt, "response")) > 0)
        xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(mf[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    weights <- model.weights(mf)
    if(!length(weights)) weights <- rep(1, nrow(x))
    fit <- rlm.default(x, y, weights, method = method,
                       wt.method = wt.method, ...)
    fit$terms <- mt
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("rlm")
    fit$call <- cl
    fit$contrasts <- attr(x, "contrasts")
    fit$xlevels <- .getXlevels(mt, mf)
    fit$na.action <- attr(mf, "na.action")
    if(model) fit$model <- mf
    if(!x.ret) fit$x <- NULL
    if(y.ret) fit$y <- y
    fit
}

rlm.default <-
  function(x, y, weights, ..., w = rep(1, nrow(x)),
           init = "ls", psi = psi.huber,
           scale.est = c("MAD", "Huber", "proposal 2"), k2 = 1.345,
           method = c("M", "MM"), wt.method = c("inv.var", "case"),
           maxit = 20, acc = 1e-4, test.vec = "resid")
{
    irls.delta <- function(old, new)
        sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
    irls.rrxwr <- function(x, w, r)
    {
        w <- sqrt(w)
        max(abs((matrix(r * w, 1, length(r)) %*% x)/
                sqrt(matrix(w, 1, length(r)) %*% (x^2))))/sqrt(sum(w * r^2))
    }

    method <- match.arg(method)
    wt.method <- match.arg(wt.method)
    nmx <- deparse(substitute(x))
    if(is.null(dim(x))) {
        x <- as.matrix(x)
        colnames(x) <- nmx
    } else x <- as.matrix(x)
    if(is.null(colnames(x)))
        colnames(x) <- paste("X", seq(ncol(x)), sep="")
    if(qr(x)$rank < ncol(x))
        stop("x is singular: singular fits are not implemented in rlm")

    if(!(any(test.vec == c("resid", "coef", "w", "NULL"))
         || is.null(test.vec))) stop("invalid testvec")
    ## deal with weights
    xx <- x
    if(!missing(weights)) {
        if(length(weights) != nrow(x))
            stop("Length of weights must equal number of observations")
        if(any(weights < 0)) stop("Negative weights value")
        if(wt.method == "inv.var") {
            fac <- sqrt(weights)
            y <- y*fac; x <- x* fac
            wt <- NULL
        } else {
            w <- w * weights
            wt <- weights
        }
    } else wt <- NULL

    if(method == "M") {
        scale.est <- match.arg(scale.est)
        if(!is.function(psi)) psi <- get(psi, mode="function")
        ## match any ... args to those of psi.
        arguments <- list(...)
        if(length(arguments)) {
            pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
            if(any(pm == 0)) warning(paste("some of ... do not match"))
            pm <- names(arguments)[pm> 0]
            formals(psi)[pm] <- unlist(arguments[pm])
        }
        if(is.character(init)) {
            temp <- if(init == "ls") lm.wfit(x, y, w, method="qr")
            else if(init == "lts") lqs(x, y, intercept=FALSE, nsamp=200)
            else stop("init method is unknown")
            coef <- temp$coef
            resid <- temp$resid
        } else {
            if(is.list(init)) coef <- init$coef
            else coef <- init
            resid <- y - x %*% coef
        }
    } else if(method == "MM") {
        scale.est <- "MM"
        temp <- lqs(x, y, intercept=FALSE, method="S", k0 = 1.548)
        coef <- temp$coef
        resid <- temp$resid
        psi <- psi.bisquare
        if(length(arguments <- list(...)))
            if(match("c", names(arguments),
                     nomatch = FALSE)) {
                c0 <- arguments$c
                if (c0 > 1.548) {
                    formals(psi)$c <- c0
                } else warning("c must be at least 1.548 and has been ignored")
            }
        scale <- temp$scale
    } else stop("method is unknown")

    done <- FALSE
    conv <- NULL
    n1 <- (if(is.null(wt)) nrow(x) else sum(wt)) - ncol(x)
    if(scale.est != "MM") scale <- mad(resid, 0)
    theta <- 2*pnorm(k2)-1
    gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
    for(iiter in 1:maxit) {
        if(!is.null(test.vec)) testpv <- get(test.vec)
        if(scale.est != "MM") {
            if(scale.est == "MAD") scale <- median(abs(resid))/0.6745
            else scale <- if(is.null(wt))
                sqrt(sum(pmin(resid^2, (k2 * scale)^2))/(n1*gamma))
            else sqrt(sum(wt*pmin(resid^2, (k2 * scale)^2))/(n1*gamma))
            if(scale == 0) {
                done <- TRUE
                break
            }
        }
        w <- psi(resid/scale)
        if(!is.null(wt)) w <- w * weights
        temp <- lm.wfit(x, y, w, method="qr")
        coef <- temp$coef
        resid <- temp$residuals
        if(!is.null(test.vec)) convi <- irls.delta(testpv, get(test.vec))
        else convi <- irls.rrxwr(x, w, resid)
        conv <- c(conv, convi)
        done <- (convi <= acc)
        if(done) break
    }
    if(!done) warning(paste("rlm failed to converge in", maxit, "steps"))
    if(!missing(weights)) {
        tmp <- (weights != 0)
        w[tmp] <- w[tmp]/weights[tmp]
    }
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl <- match.call()
    cl[[1]] <- as.name("rlm")
    fit <- list(coefficients = coef, residuals = resid, effects = temp$effects,
                rank = temp$rank, fitted.values = temp$fitted.values,
                assign = temp$assign,  qr = temp$qr, df.residual = NA, w = w,
                s = scale, psi = psi, k2 = k2,
                weights = if(!missing(weights)) weights,
                conv = conv, converged = done, x = xx, call = cl)
    class(fit) <- c("rlm", "lm")
    fit
}

print.rlm <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    if(x$converged)
        cat("Converged in", length(x$conv), "iterations\n")
    else cat("Ran", length(x$conv), "iterations without convergence\n")
    coef <- x$coef
    cat("\nCoefficients:\n")
    print(coef, ...)
    nobs <- length(x$resid)
    rdf <- nobs - length(coef)
    cat("\nDegrees of freedom:", nobs, "total;", rdf, "residual\n")
    cat("Scale estimate:", format(signif(x$s,3)), "\n")
    invisible(x)
}

summary.rlm <- function(object, method=c("XtX", "XtWX"),
                        correlation = TRUE, ...)
{
    method <- match.arg(method)
    s <- object$s
    coef <- object$coef
    ptotal <- length(coef)
    resid <- object$resid
    n <- length(resid)
    if(any(na <- is.na(coef))) coef <- coef[!na]
    cnames <- names(coef)
    p <- length(coef)
    rdf <- n - p
    rinv <- diag(p)
    dimnames(rinv) <- list(cnames, cnames)
    w <- object$psi(resid/s)
    S <- sum((resid*w)^2)/rdf
    psiprime <- object$psi(resid/s, deriv=1)
    mn <- mean(psiprime)
    kappa <- 1 + p*var(psiprime)/(n*mn^2)
    stddev <- sqrt(S)*(kappa/mn)
    X <- object$x * sqrt(object$weights)
    if(method == "XtWX")  X <- X * sqrt(w/mean(w))
    R <- qr(X)$qr
    R <- R[1:p, 1:p, drop = FALSE]
    R[lower.tri(R)] <- 0
    rinv <- solve(R, rinv)
    dimnames(rinv) <- list(cnames, cnames)
    rowlen <- (rinv^2 %*% rep(1, p))^0.5
    names(rowlen) <- cnames
    if(correlation) {
        correl <- rinv * array(1/rowlen, c(p, p))
        correl <- correl %*% t(correl)
    } else correl <- NULL
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
    class(object) <- "summary.rlm"
    object
}

print.summary.rlm <-
function(x, digits = max(3, .Options$digits - 3), ...)
{
    cat("\nCall: ")
    dput(x$call)
    resid <- x$residuals
    df <- x$df
    rdf <- df[2]
    if(rdf > 5) {
        cat("Residuals:\n")
        if(length(dim(resid)) == 2) {
            rq <- apply(t(resid), 1, quantile)
            dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q", "Max"),
                                 colnames(resid))
        } else {
            rq <- quantile(resid)
            names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
        }
        print(rq, digits = digits, ...)
    } else if(rdf > 0) {
        cat("Residuals:\n")
        print(resid, digits = digits, ...)
    }
    if(nsingular <- df[3] - df[1])
        cat("\nCoefficients: (", nsingular,
            " not defined because of singularities)\n", sep = "")
    else cat("\nCoefficients:\n")
    print(format(round(x$coef, digits = digits)), quote = FALSE, ...)
    cat("\nResidual standard error:", format(signif(x$sigma, digits)),
        "on", rdf, "degrees of freedom\n")
    if(!is.null(correl <- x$correlation)) {
        p <- dim(correl)[2]
        if(p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1, -p, drop = FALSE], quote = FALSE, digits = digits, ...)
        }
    }
    invisible(x)
}

psi.huber <- function(u, k = 1.345, deriv=0)
{
    if(!deriv) return(pmin(1, k / abs(u)))
    abs(u) <= k
}

psi.hampel <- function(u, a = 2, b = 4, c = 8, deriv=0)
{
    U <- pmin(abs(u) + 1e-50, c)
    if(!deriv) return(ifelse(U <= a, U, ifelse(U <= b, a, a*(c-U)/(c-b) ))/U)
    ifelse(abs(u) <= c, ifelse(U <= a, 1, ifelse(U <= b, 0, -a/(c-b))), 0)
}

psi.bisquare <- function(u, c = 4.685, deriv=0)
{
    if(!deriv) return((1  - pmin(1, abs(u/c))^2)^2)
    t <- (u/c)^2
    ifelse(t < 1, (1 - t)*(1 - 5*t), 0)
}

se.contrast.rlm <-
    function(object, contrast.obj, coef = contr.helmert(ncol(contrast))[, 1],
             data = NULL, ...)
{
    contrast.weight.aov <- function(object, contrast)
    {
        asgn <- object$assign[object$qr$pivot[1:object$rank]]
        uasgn <- unique(asgn)
        nterms <- length(uasgn)
        nmeffect <- c("(Intercept)",
                      attr(object$terms, "term.labels"))[1 + uasgn]
        effects <- as.matrix(qr.qty(object$qr, contrast))
        res <- matrix(0, nrow = nterms, ncol = ncol(effects),
                      dimnames = list(nmeffect, colnames(contrast)))
        for(i in seq(nterms)) {
            select <- (asgn == uasgn[i])
            res[i,] <- colSums(effects[seq(along=asgn)[select], , drop = FALSE]^2)
        }
        res
    }
    if(is.null(data)) contrast.obj <- eval(contrast.obj)
    else contrast.obj <- eval(substitute(contrast.obj), data, parent.frame())
    if(!is.matrix(contrast.obj)) { # so a list
        if(sum(coef) != 0)
            stop("coef must define a contrast, i.e., sum to 0")
        if(length(coef) != length(contrast.obj))
            stop("coef must have same length as contrast.obj")
        contrast <-
            sapply(contrast.obj, function(x)
               {
                   if(!is.logical(x))
                       stop(paste("Each element of", substitute(contrasts.list),
                                  " must be\nlogical"))
                   x/sum(x)
               })
        contrast <- contrast %*% coef
        if(!any(contrast) || all(is.na(contrast)))
            stop("The contrast defined is empty (has no TRUE elements)")
    } else {
        contrast <- contrast.obj
        if(any(abs(colSums(contrast)) > 1e-8))
            stop("Columns of contrast.obj must define a contrast (sum to zero)")
        if(length(colnames(contrast)) == 0)
            colnames(contrast) <- paste("Contrast", seq(ncol(contrast)))
    }
    weights <- contrast.weight.aov(object, contrast)
    object$sigma * if(!is.matrix(contrast.obj)) sqrt(sum(weights)) else sqrt(colSums(weights))
}

predict.rlm <- function (object, newdata = NULL, scale = NULL, ...)
{
    ## problems with using predict.lm are the scale and
    ## the QR decomp which has been done on down-weighted values.
    object$qr <- qr(sqrt(object$weights) * object$x)
    predict.lm(object, newdata = newdata, scale = object$s, ...)
}
