spec.taper <- function(x, p = 0.1)
{
    if (any(p < 0) || any(p > 0.5))
        stop("p must be between 0 and 0.5")
    x <- as.ts(x)
    a <- attributes(x)
    x <- as.matrix(x)
    nc <- ncol(x)
    if (length(p) == 1)
        p <- rep(p, nc)
    else if (length(p) != nc)
        stop("length of p must be 1 or equal the number of columns of x")
    nr <- nrow(x)
    for (i in 1 : nc) {
        m <- floor(nr * p[i])
        w <- 0.5 * (1 - cos(pi * seq(1, 2 * m - 1, by = 2) / (2 * m)))
        x[, i] <- c(w, rep(1, nr - 2 * m), rev(w)) * x[, i]
    }
    attributes(x) <- a
    x
}

wt.huber <- function(u, c = 1.345)
    ifelse(abs(u) < c, 1, c / abs(u))

.First.lib <- function(lib, pkg)
{
    if(version$major==0 && version$minor < 63)
        stop("This version for R 0.63 or later")
    library.dynam("MASS", pkg, lib)
}

profile <- function(fitted, ...) UseMethod("profile")
