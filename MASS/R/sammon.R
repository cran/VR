# file MASS/sammon.q
# copyright (C) 1994-2003 W. N. Venables and B. D. Ripley
#
sammon <- function(d, y= cmdscale(d, k), k=2, niter=100, trace=TRUE,
                   magic=0.2, tol=1e-4)
{
    call <- match.call()
    if(any(is.infinite(as.vector(d))))
        stop("Infs not allowed in 'd'")
    if(any(is.na(d)) && missing(y))
        stop("an initial configuration must be supplied if there are NAs in 'd'")
    if(is.null(n <- attr(d, "Size"))) {
        x <- as.matrix(d)
        if((n <- nrow(x)) != ncol(x))
            stop("distances must be result of dist or a square matrix")
    }
    else {
        x <- matrix(0, n, n)
        x[row(x) > col(x)] <- d
        x <- x + t(x)
    }
    ab <- x[row(x) < col(x)]<=0
    if (any(ab, na.rm = TRUE)) {
        ab <- !is.na(ab) & ab
        aa <- cbind(as.vector(row(x)), as.vector(col(x)))[row(x) < col(x),]
        aa <- aa[ab, , drop=FALSE]
        stop(sprintf(gettext("zero or negative distance between objects %d and %d"),
                     aa[1,1], aa[1,2]), domain = NA)
    }
    nas <- is.na(x)
    diag(nas) <- FALSE  # diag never used
    if(any(rowSums(!nas) < 2)) stop("not enough non-missing data")

    if(!is.matrix(y)) stop("'y' must be a matrix")
    if(any(dim(y) != c(n, k)) ) stop("invalid initial configuration")
    if(any(!is.finite(y))) stop("initial configuration must be complete")
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    z <- .C("VR_sammon",
            x = x,
            as.integer(n),
            as.integer(k),
            y = y,
            as.integer(niter),
            e = double(1),
            as.integer(trace),
            as.double(magic),
            as.double(tol),
            NAOK = TRUE, PACKAGE = "MASS"
            )
    points <- z$y
    rn <- if(is.matrix(d)) rownames(d) else names(d)
    dimnames(points) <- list(rn, NULL)
    list(points=points, stress=z$e, call=call)
}
