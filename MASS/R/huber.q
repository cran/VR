# file MASS/huber.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
huber <- function(y, k=1.5, tol = 1.0e-6)
{
    y <- y[!is.na(y)]
    n <- length(y)
    mu <- median(y)
    s <- mad(y)
    repeat{
        yy <- pmin(pmax(mu-k*s,y),mu+k*s)
        mu1 <- sum(yy)/n
        if(abs(mu-mu1) < tol*s) break
        mu <- mu1
    }
    list(mu=mu,s=s)
}
