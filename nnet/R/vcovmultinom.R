# vcov.multinom <- function(object, ...)
# {
#     hess <- multinomHess(object)
#     structure(chol2inv(chol(hess)), dimnames = dimnames(hess))
# }

multinomHess <- function(object)
{
    probs <- fitted(object)
    Z <- model.matrix(object)
    coefs <- coef(object)
    if (is.vector(coefs)){ # ie there are only 2 response categories
        coefs <- t(as.matrix(coefs))
        probs <- cbind(1 - probs, probs)
    }
    coefdim <- dim(coefs)
    p <- coefdim[2]
    k <- coefdim[1]
    ncoefs <- k * p
    kpees <- rep(p, k)
    n <- dim(Z)[1]
##  Now compute the observed (= expected, in this case) information,
##  e.g. as in T Amemiya "Advanced Econometrics" (1985) pp295-6.
##  Here i and j are as in Amemiya, and x, xbar are vectors
##  specific to (i,j) and to i respectively.
    info <- matrix(0, ncoefs, ncoefs)
    Names <- dimnames(coefs)
    if (is.null(Names[[1]])) Names <- Names[[2]]
    else Names <- as.vector(outer(Names[[2]], Names[[1]],
                             function(name2, name1)
                                 paste(name1, name2, sep = ":")))
    dimnames(info) <- list(Names, Names)
    x0 <- matrix(0, p, k+1)
    row.totals <- object$weights
    for (i in 1:n){
        Zi <- Z[i, ]
        xbar <- Zi * rep(probs[i, -1, drop=FALSE], kpees)
        for (j in 1:(k+1)){
            x <- x0
            x[, j] <- Zi
            x <- x[, -1, drop = FALSE]
            x <- x - xbar
            dim(x) <- c(1, ncoefs)
            info <- info + (row.totals[i] * probs[i, j] * crossprod(x))
        }
    }
    info
}
