# file MASS/misc.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
Choleski <- function(x)
{
    .NotYetImplemented()
    x <- as.Matrix(x)
    class(x) <- Matrix.class(x)
    if(!inherits(x, "Hermitian")) stop("x must be symmetric")
    LU <- expand(lu(x))
    B <- LU$block.diagonal
    if(!inherits(B, "Diagonal") || any(diag(B) <= 0)
       || !inherits(LU$permutation, "Identity"))
        stop("x is not positive definite")
    x <- LU$triangular %*% sqrt(B)
    class(x) <- Matrix.class(x)
    x
}

mat2tr <- function(z)
{
    dn <- names(dimnames(z))
    dx <- dimnames(z)[[1]]
    x <- as.numeric(substring(dx, nchar(dn[1]) + 2))
    dy <- dimnames(z)[[2]]
    y <- as.numeric(substring(dy, nchar(dn[2]) + 2))
    cbind(expand.grid(x = x, y = y), z = as.vector(z))
}

con2tr <- function(obj)
{
    data.frame(expand.grid(x=obj$x,y=obj$y),z=as.vector(obj$z))
}

Null <- function(M)
{
    tmp <- qr(M)
    set <- if(tmp$rank == 0) 1:ncol(M) else  - (1:tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}

ginv <- function(X, tol = sqrt(.Machine$double.eps))
{
    #
    # based on suggestions of R M Heiberger, T M Hesterburg and WNV
    #
    if(length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if(!is.matrix(X)) X <- as.matrix(X)
    Xsvd <- svd(X)
    if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if(!any(Positive)) array(0, dim(X)[2:1])
    else Xsvd$v[, Positive] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive]))
}
