# file MASS/hubers.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
hubers <- function(y, k=1.5, mu, s, initmu=median(y), tol = 1.0e-6)
{
    y <- y[!is.na(y)]
    n <- length(y)
    if(missing(mu)) {
        mu0 <- initmu
        n1 <- n - 1
    } else {
        mu0 <- mu
        mu1 <- mu
        n1 <- n
    }
    if(missing(s)) s0 <- mad(y) else {s0 <- s;s1 <- s}
    th <- 2*pnorm(k)-1
    beta <- th + k^2*(1-th) - 2*k*dnorm(k)
    repeat{
        yy <- pmin(pmax(mu0-k*s0,y), mu0+k*s0)
        if(missing(mu)) mu1 <- sum(yy)/n
        if(missing(s)) {
            ss <- sum((yy-mu1)^2)/n1
            s1 <- sqrt(ss/beta)
        }
        if((abs(mu0-mu1) < tol*s0) && abs(s0-s1) < tol*s0) break
        mu0 <- mu1; s0 <- s1
    }
    list(mu=mu0, s=s0)
}
