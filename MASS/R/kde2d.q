# file MASS/kde2d.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
kde2d <- function(x, y, h, n = 25, lims=c(range(x), range(y)) )
{
    nx <- length(x)
    if(length(y) != nx)
        stop("Data vectors must be the same length")
    gx <- seq(lims[1], lims[2], length = n)
    gy <- seq(lims[3], lims[4], length = n)
    if (missing(h))
        h <- c(bandwidth.nrd(x), bandwidth.nrd(y))
    h <- h/4                            # for S's bandwidth scale
    ax <- outer(gx, x, "-" )/h[1]
    ay <- outer(gy, y, "-" )/h[2]
    z <- matrix(dnorm(ax), n, nx) %*%
        t(matrix(dnorm(ay),n, nx))/ (nx * h[1] * h[2])
    return(list(x = gx, y = gy, z = z))
}
