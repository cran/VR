# file MASS/hist.scott.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#

hist.scott <- function(x, prob = TRUE, xlab = deparse(substitute(x)), ...)
   invisible(hist(x, nclass.scott(x), prob=prob, xlab=xlab, ...))
hist.FD <- function(x, prob = TRUE, xlab = deparse(substitute(x)), ...)
   invisible(hist(x, nclass.FD(x), prob=prob, xlab=xlab, ...))

frequency.polygon <- function(x, nclass = nclass.freq(x),
    xlab="", ylab="", ...)
{
    hst <- hist(x, nclass, probability=TRUE, plot=FALSE, ...)
    midpoints <- 0.5 * (hst$breaks[-length(hst$breaks)]
                        + hst$breaks[-1])
    plot(midpoints, hst$counts, type="l", xlab=xlab, ylab=ylab)
}

nclass.freq <- function(x)
{
    h <- 2.15 * sqrt(var(x)) * length(x)^(-1/5)
    ceiling(diff(range(x))/h)
}

bandwidth.nrd <- function(x)
{
    r <- quantile(x, c(0.25, 0.75))
    h <- (r[2] - r[1])/1.34
    4 * 1.06 * min(sqrt(var(x)), h) * length(x) ^ (-1/5)
}
