# file MASS/lm.ridge.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
lm.ridge <- function(formula, data, subset, na.action,
    lambda = 0, model = FALSE, x = FALSE, y = FALSE, contrasts = NULL, ...)
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$model <- m$x <- m$y <- m$contrasts <- m$... <- m$lambda <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    Y <- model.extract(m, response)
    X <- model.matrix(Terms, m, contrasts)
    n <- nrow(X); p <- ncol(X)
    if(Inter <- attr(Terms, "intercept"))
    {
        Xm <- apply(X[, -Inter], 2, mean)
        Ym <- mean(Y)
        p <- p - 1
        X <- X[, -Inter] - rep(Xm, rep(n, p))
        Y <- Y - Ym
    } else Ym <- Xm <- NA
    Xscale <- drop(rep(1/n, n) %*% X^2)^0.5
    X <- X/rep(Xscale, rep(n, p))
    Xs <- svd(X)
    rhs <- t(Xs$u) %*% Y
    d <- Xs$d
    lscoef <-  Xs$v %*% (rhs/d)
    lsfit <- X %*% lscoef
    resid <- Y - lsfit
    s2 <- sum(resid^2)/(n - p - Inter)
    HKB <- (p-2)*s2/sum(lscoef^2)
    LW <- (p-2)*s2*n/sum(lsfit^2)
    k <- length(lambda)
    div <- d^2 + rep(lambda, rep(p,k))
    a <- drop(d*rhs)/div
    dim(a) <- c(p, k)
    coef <- Xs$v %*% a
    dimnames(coef) <- list(names(Xscale), format(lambda))
    GCV <- apply((Y - X %*% coef)^2, 2, sum)/
        (n-apply(matrix(d^2/div,p), 2, sum))^2
    res <- list(coef = drop(coef), scales = Xscale,
                Inter = Inter, lambda = lambda, ym = Ym, xm = Xm,
                GCV = GCV, kHKB = HKB, kLW = LW)
    class(res) <- "ridgelm"
    res
}

print.ridgelm <- function(x, ...)
{
    scaledcoef <- t(as.matrix(x$coef / x$scales))
    if(x$Inter) {
        inter <- x$ym - scaledcoef %*% x$xm
        scaledcoef<- cbind(Intercept=inter, scaledcoef)
    }
    print(drop(scaledcoef), ...)
}

select <- function(obj) UseMethod("select")

select.ridgelm <- function(obj)
{
    cat("modified HKB estimator is", format(obj$kHKB), "\n")
    cat("modified L-W estimator is", format(obj$kLW), "\n")
    GCV <- obj$GCV
    if(length(GCV) > 0) {
        k <- seq(along=GCV)[GCV==min(GCV)]
        cat("smallest value of GCV  at",
            format(obj$lambda[k]), "\n")
    }
}

plot.ridgelm <- function(x)
    matplot(x$lambda, t(x$coef), type = "l")
