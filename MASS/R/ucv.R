# file MASS/ucv.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#

width.SJ <- function(x, nb=1000, lower=0.1*hmax, upper=hmax,
		     method = c("ste", "dpi"))
{
    fSD <- function(h, x, alph2, c1, n, d)
        (c1/SDh(x, alph2 * h^(5/7), n, d))^(1/5) - h
    SDh <- function(x, h, n, d)
        .C("VR_phi4_bin",
           as.integer(n),
           as.integer(length(x)),
           as.double(d),
           x,
           as.double(h),
           u = double(1))$u
    TDh <- function(x, h, n, d)
        .C("VR_phi6_bin",
           as.integer(n),
           as.integer(length(x)),
           as.double(d),
           x,
           as.double(h),
           u = double(1))$u

    method <- match.arg(method)
    n <- length(x)
    if(!is.loaded(symbol.C("VR_phi4_bin")))
        stop("Compiled code has not been dynamically loaded")
    storage.mode(x) <- "double"
    n <- length(x)
    Z <- .C("VR_den_bin",
            as.integer(n),
            as.integer(nb),
            d = double(1),
            x,
            cnt = integer(nb)
            )
    d <- Z$d; cnt <- as.integer(Z$cnt)
    hmax <- 1.144 * sqrt(var(x)) * n^(-1/5)
    scale <- min(sqrt(var(x)), IQR(x)/1.349)
    a <- 1.24 * scale * n^(-1/7)
    b <- 1.23 * scale * n^(-1/9)
    c1 <- 1/(2*sqrt(pi)*n)
    TD  <- -TDh(cnt, b, n, d)
    alph2 <- 1.357*(SDh(cnt, a, n, d)/TD)^(1/7)
    if(method == "dpi")
        res <- (c1/SDh(cnt,(2.394/(n * TD))^(1/7) , n, d))^(1/5)
    else {
        if (fSD(lower, cnt, alph2, c1, n, d) *
            fSD(upper, cnt, alph2, c1, n, d) > 0)
            stop("No solution in the specified range of bandwidths")
        res <- uniroot(fSD, c(lower, upper), tol=0.1*lower,
                       x=cnt, alph2=alph2, c1=c1, n=n, d=d)$root
    }
    4 * res
}


ucv <- function(x, nb=1000, lower=0.1*hmax, upper=hmax)
{
    fucv <- function(h, x, n, d)
        .C("VR_ucv_bin",
           as.integer(n),
           as.integer(length(x)),
           as.double(d),
           x,
           as.double(h),
           u = double(1))$u

    n <- length(x)
    hmax <- 1.144 * sqrt(var(x)) * n^(-1/5) * 4
    if(!is.loaded(symbol.C("VR_den_bin")))
        stop("Compiled code has not been dynamically loaded")
    storage.mode(x) <- "double"
    Z <- .C("VR_den_bin",
            as.integer(n),
            as.integer(nb),
            d = double(1),
            x,
            cnt = integer(nb)
            )
    d <- Z$d; cnt <- as.integer(Z$cnt)
    h <- optimize(fucv, c(lower, upper), tol=0.1*lower,
                  x=cnt, n=n, d=d)$minimum
    if(h < 1.1*lower | h > upper-0.1*lower)
        warning("minimum occurred at one end of the range")
    h
}


bcv <- function(x, nb=1000, lower=0.1*hmax, upper=hmax)
{
    fbcv <- function(h, x, n, d)
        .C("VR_bcv_bin",
           as.integer(n),
           as.integer(length(x)),
           as.double(d),
           x,
           as.double(h),
           u = double(1))$u

    n <- length(x)
    hmax <- 1.144 * sqrt(var(x)) * n^(-1/5) * 4
    if(!is.loaded(symbol.C("VR_den_bin")))
        stop("Compiled code has not been dynamically loaded")
    storage.mode(x) <- "double"
    Z <- .C("VR_den_bin",
            as.integer(n),
            as.integer(nb),
            d = double(1),
            x,
            cnt = integer(nb)
            )
    d <- Z$d; cnt <- as.integer(Z$cnt)
    h<- optimize(fbcv, c(lower, upper), tol=0.1*lower,
                 x=cnt, n=n, d=d)$minimum
    if(h < 1.1*lower | h > upper-0.1*lower)
        warning("minimum occurred at one end of the range")
    h
}
