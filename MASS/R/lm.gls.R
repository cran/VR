# file MASS/lm.gls.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
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
    Y <- model.extract(m, response)
    X <- model.matrix(Terms, m, contrasts)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    n <- nrow(X)
    if(any(dim(W) != c(n, n))) stop("dim(W) is not correct")
    eW <- eigen(W, TRUE)
    d <- eW$values
    if(any(d <= 0)) stop("W is not positive definite")
    A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
    fit <- lm.fit(A %*% X, A %*% Y, method=method, ...)
    fit$terms <- Terms
    fit$call <- call
    if(model) fit$model <- m
    if(x) fit$x <- X
    if(y) fit$y <- Y
    attr(fit, "na.message") <- attr(m, "na.message")
    if(!is.null(attr(m, "na.action"))) fit$na.action <- attr(m, "na.action")
    class(fit) <- c("lm.gls", class(fit))
    fit$xlevels <- xlev
    fit$contrasts <- attr(X, "contrasts")
    fit
}
