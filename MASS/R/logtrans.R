# file MASS/logtrans.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
logtrans <- function(object, ...) UseMethod("logtrans")

logtrans.default<-
function(object, ..., alpha = seq(0.5, 6, by = 0.25) - min(y),
	plotit = TRUE, interp = (plotit && (m <
	100)), xlab = "alpha", ylab = "log Likelihood")
{
    if(is.null(object$y) || is.null(object$qr))
        stop(paste(deparse(substitute(object)),
                   "does not have both 'qr' and 'y' components"))
    y <- object$y
    n <- length(y)
    if(any(y + min(alpha) <= 0))
        stop("Response variable must be positive after additions")
    xqr <- object$qr
    xl <- loglik <- as.vector(alpha)
    m <- length(xl)
    for(i in 1:m) {
        rs <- qr.resid(xqr, yt <- log(y + alpha[i]))
        loglik[i] <-  - n/2 * log(sum(rs^2)) - sum(yt)
    }
    if(interp) {
        sp <- spline(alpha, loglik, n = 100)
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
    invisible(list(x = xl, y = loglik))
}

logtrans.formula <-
function(object, data = NULL, ...)
{
  object <- aov(object, data = data, y = TRUE, qr = TRUE)
  invisible(NextMethod("logtrans"))
}

logtrans.lm <- function(object, ...)
{
    if(is.null(object$y) || is.null(object$qr))
        object <- update(object, y = TRUE, qr = TRUE)
    invisible(NextMethod("logtrans"))
}
