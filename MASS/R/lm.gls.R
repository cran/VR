# file MASS/lm.gls.q
# copyright (C) 1994-2005 W. N. Venables and B. D. Ripley
#
lm.gls <-
    function(formula, data, W, subset, na.action, inverse = FALSE,
             method = "qr",
             model = FALSE, x = FALSE, y = FALSE, contrasts = NULL, ...)
{
    call <- match.call()
    m <- match.call(expand = FALSE)
    m$W <- m$inverse <- m$method <- m$model <- m$x <-
        m$y <- m$contrasts <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    if(method == "model.frame") return(m)
    Terms <- attr(m, "terms")
    Y <- model.response(m)
    X <- model.matrix(Terms, m, contrasts)
    n <- nrow(X)
    if(any(dim(W) != c(n, n))) stop("dim(W) is not correct")
    eW <- eigen(W, TRUE)
    d <- eW$values
    if(any(d <= 0)) stop("'W' is not positive definite")
    A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
    Ainv <- eW$vector %*% diag(d^ifelse(inverse, 0.5, -0.5))
    fit <- lm.fit(A %*% X, A %*% Y, method=method, ...)
    ## revert to original coords
    fit$fitted.values <- drop(Ainv %*% fit$fitted.values)
    fit$residuals <- drop(Ainv %*% fit$residuals)
    fit$terms <- Terms
    fit$call <- call
    if(model) fit$model <- m
    if(x) fit$x <- X
    if(y) fit$y <- Y
    fit$na.action <- attr(m, "na.action")
    class(fit) <- "lm.gls"
    fit$xlevels <- .getXlevels(Terms, m)
    fit$contrasts <- attr(X, "contrasts")
    fit
}
