# file MASS/negexp.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
negexp.ival <- function(x, b0, b1, th)
{
    pnames <- as.character(sys.call()[3:5])
    y <- get(".nls.initial.response")
    if(length(unique(x)) < 3)
        stop("at least 3 distinct x values are needed")
    mx <- mean(x)
    b <- as.vector(lsfit(cbind(x - mx, -(x - mx)^2/2), y)$coef)
    rx <- range(x)
    xh <- mx + b[2]/b[3]
    if(prod(xh - rx) < 0)
        if(xh - rx[1] > rx[2] - xh) rx[2] <- xh  else rx[1] <- xh
    x0 <- c(rx[1], sum(rx)/2, rx[2])
    dy <- diff(b[1] + b[2]*(x0 - mx) - (b[3]*(x0 - mx)^2)/2)
    th <- (x0[2] - x0[1])/log(dy[1]/dy[2])
    b <- as.vector(lsfit(exp( - x/th), y)$coef)
    pars <- list(b[1], b[2], th)
    names(pars) <- pnames
    print(unlist(pars))
    pars
}

negexp.SSival <- function(mCall, data, LHS)
{
    x <- eval(mCall[["x"]], data)
    if(length(x) < 3)
        stop("at least 3 distinct x values are needed")
    y <- eval(LHS, data)
    mx <- mean(x)
    b <- as.vector(lsfit(cbind(x - mx,  - (x - mx)^2/2), y)$coef)
    rx <- range(x)
    xh <- mx + b[2]/b[3]
    if(prod(xh - rx) < 0)
        if(xh - rx[1] > rx[2] - xh)
            rx[2] <- xh
        else rx[1] <- xh
    x0 <- c(rx[1], sum(rx)/2, rx[2])
    dy <- diff(b[1] + b[2] * (x0 - mx) - (b[3] * (x0 - mx)^2)/2)
    th <- (x0[2] - x0[1])/log(dy[1]/dy[2])
    b <- as.vector(lsfit(exp( - x/th), y)$coef)
    pars <- list(b[1], b[2], th)
    names(pars) <- mCall[c("b0", "b1", "th")]
    print(unlist(pars))
    pars
}
