# file MASS/gamma.shape.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
gamma.shape <- function(object, ...) UseMethod("gamma.shape")

gamma.shape.glm <- function(object, it.lim = 10,
			    eps.max = .Machine$double.eps^0.25,
			    verbose = FALSE)
{
    if(is.null(object$y)) object <- update(object, y = TRUE)
    y <- object$y
    A <- object$prior.weights
    if(is.null(A)) A <- rep(1, length(y))
    u <- fitted(object)
    Dbar <- object$deviance/object$df.residual
    alpha <- (6 + 2*Dbar)/(Dbar*(6 + Dbar))
    if(verbose) cat("Initial estimate:", format(alpha), "\n")
    fixed <-  -y/u - log(u) + log(A) + 1 + log(y + (y == 0))
    eps <- 1
    itr <- 0
    while(abs(eps) > eps.max && (itr <- itr + 1) <= it.lim) {
        sc <- sum(A * (fixed + log(alpha) - digamma(A * alpha)))
        inf <- sum(A * (A * trigamma(A * alpha) - 1/alpha))
        alpha <- alpha + (eps <- sc/inf)
        if(verbose) cat("Iter. ", itr, " Alpha:", alpha, "\n")
    }
    if(itr > it.lim) warning("Iteration limit reached")
    res <- list(alpha = alpha, SE = sqrt(1/inf))
    class(res) <- "gamma.shape"
    res
}

gamma.dispersion <- function(object, ...)
    1/gamma.shape(object, ...)[[1]]

print.gamma.shape <- function(x, ...)
{
    y <- x
    x <- array(unlist(x), dim = 2:1,
               dimnames = list(c("Alpha:", "SE:"), ""))
    NextMethod("print")
    invisible(y)
}
