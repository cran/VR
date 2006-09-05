# file MASS/boxcox.q
# copyright (C) 1994-2004 W. N. Venables and B. D. Ripley
#
boxcox <- function(object, ...) UseMethod("boxcox")

boxcox.formula <-
function(object, lambda = seq(-2, 2, 1/10), plotit =  TRUE,
         interp = (plotit && (m < 100)), eps = 1/50,
         xlab = expression(lambda), ylab = "log-Likelihood", ...)
{
    m <- length(lambda)
    object <- lm(object, y = TRUE, qr = TRUE, ...)
    result <- NextMethod()
    if(plotit) invisible(result)
    else result
}

boxcox.lm <-
function(object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
         interp = (plotit && (m < 100)), eps = 1/50,
         xlab = expression(lambda), ylab = "log-Likelihood", ...)
{
    if(is.null(object$y) || is.null(object$qr))
        object <- update(object, y = TRUE, qr = TRUE, ...)
    result <- NextMethod()
    if(plotit) invisible(result)
    else result
}

boxcox.default <-
function(object, lambda = seq(-2, 2, 1/10), plotit = TRUE,
         interp = (plotit && (m < 100)), eps = 1/
         50, xlab = expression(lambda), ylab = "log-Likelihood", ...)
{
    if(is.null(object$y) || is.null(object$qr))
        stop(paste(deparse(substitute(object)),
                   "does not have both 'qr' and 'y' components"
                   ))
    y <- object$y
    n <- length(y)
    if(any(y <= 0))
        stop("response variable must be positive")
    xqr <- object$qr
    logy <- log(y)
    ydot <- exp(mean(logy))
    xl <- loglik <- as.vector(lambda)
    m <- length(xl)
    for(i in 1:m) {
        if(abs(la <- xl[i]) > eps)
            yt <- (y^la - 1)/la
        else yt <- logy * (1 + (la * logy)/2 *
                           (1 + (la * logy)/3 * (1 + (la * logy)/4)))
        loglik[i] <-  - n/2 * log(sum(qr.resid(xqr, yt/ydot^(la - 1))^2))
    }
    if(interp) {
        sp <- spline(xl, loglik, n = 100)
        xl <- sp$x
        loglik <- sp$y
        m <- length(xl)
    }
    if(plotit) {
        mx <- (1:m)[loglik == max(loglik)][1]
        Lmax <- loglik[mx]
        lim <- Lmax - qchisq(19/20, 1)/2
        plot(xl, loglik, xlab = xlab, ylab = ylab, type
             = "l", ylim = range(loglik, lim))
        plims <- par("usr")
        abline(h = lim, lty = 3)
        y0 <- plims[3]
        scal <- (1/10 * (plims[4] - y0))/par("pin")[2]
        scx <- (1/10 * (plims[2] - plims[1]))/par("pin")[1]
        text(xl[1] + scx, lim + scal, " 95%")
        la <- xl[mx]
        if(mx > 1 && mx < m)
            segments(la, y0, la, Lmax, lty = 3)
        ind <- range((1:m)[loglik > lim])
        if(loglik[1] < lim) {
            i <- ind[1]
            x <- xl[i - 1] + ((lim - loglik[i - 1]) *
                              (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
            segments(x, y0, x, lim, lty = 3)
        }
        if(loglik[m] < lim) {
            i <- ind[2] + 1
            x <- xl[i - 1] + ((lim - loglik[i - 1]) *
                              (xl[i] - xl[i - 1]))/(loglik[i] - loglik[i - 1])
            segments(x, y0, x, lim, lty = 3)
        }
    }
    list(x = xl, y = loglik)
}
