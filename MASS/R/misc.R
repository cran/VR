# file MASS/misc.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#

Null <- function(M)
{
    tmp <- qr(M)
    set <- if(tmp$rank == 0) 1:ncol(M) else  - (1:tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}

ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
#
# based on suggestions of R. M. Heiberger, T. M. Hesterberg and WNV
#
    if(length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if(!is.matrix(X)) X <- as.matrix(X)
    Xsvd <- svd(X)
    if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if(!any(Positive)) array(0, dim(X)[2:1])
    else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
}
