# file MASS/neg.bin.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
neg.bin <- function (theta = stop("theta must be given"))
{
    .Theta <- theta
    stats <- make.link("log")
    variance <- function(mu)
        mu + mu^2/.Theta
    validmu <- function(mu)
        all(mu > 0)
    dev.resids <- function(y, mu, wt)
        2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) *
                  log((y + .Theta)/ (mu + .Theta)))
    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Theta) * log((y + .Theta)/ (mu + .Theta)) - y * log(mu) +
            lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) - lgamma(.Theta+y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop(paste("Negative values not allowed for",
                                   "the Poisson family"))
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
    structure(list(family = "Negative Binomial", link = "log",
                   linkfun = stats$linkfun, linkinv = stats$linkinv,
                   variance = variance, dev.resids = dev.resids,
                   aic = aic, mu.eta = stats$mu.eta,
                   initialize = initialize, validmu = validmu,
                   valideta = stats$valideta), class = "family")
}
