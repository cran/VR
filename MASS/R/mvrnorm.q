# file MASS/mvrnorm.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
mvrnorm <- function(n = 1, mu, Sigma, tol=1e-6)
{
    p <- length(mu)
    if(!all(dim(Sigma) == c(p,p))) stop("incompatible arguments")
    eS <- eigen(Sigma, sym = TRUE)
    ev <- eS$values
    if(!all(ev >= -tol*abs(ev[1]))) stop("Sigma is not positive definite")
    X <- mu + eS$vectors %*% diag(sqrt(pmax(ev, 0))) %*% matrix(rnorm(p*n),p)
    nm <- names(mu)
    if(is.null(nm) && !is.null(dn <- dimnames(Sigma))) nm <- dn[[1]]
    dimnames(X) <- list(nm, NULL)
    if(n == 1) drop(X) else t(X)
}
