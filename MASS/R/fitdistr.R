fitdistr <- function(x, densfun, start, ...)
{
    myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
    mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))
    mydt <- function(x, m, s, df, log) dt((x-m)/s, df, log = TRUE) - log(s)

    if(missing(start)) start <- NULL
    dots <- names(list(...))
    dots <- dots[!is.element(dots, c("upper", "lower"))]
    if(missing(x) || length(x) == 0 || mode(x) != "numeric")
        stop("'x' must be a non-empty numeric vector")
    if(missing(densfun) || !(is.function(densfun) || is.character(densfun)))
        stop("'densfun' must be supplied as a function or name")
    if(is.character(densfun)) {
        distname <- tolower(densfun)
        densfun <-
            switch(distname,
                   "beta" = dbeta,
                   "cauchy" = dcauchy,
                   "chi-squared" = dchisq,
                   "exponential" = dexp,
                   "f" = df,
                   "gamma" = dgamma,
                   "log-normal" = dlnorm,
                   "lognormal" = dlnorm,
                   "logistic" = dlogis,
                   "negative binomial" = dnbinom,
                   "normal" = dnorm,
                   "t" = mydt,
                   "uniform" = dunif,
                   "weibull" = dweibull,
                   NULL)
        if(is.null(densfun)) stop("unsupported distribution")
        if(distname == "normal") {
            if(!is.null(start))
                stop("supplying pars for the Normal is not supported")
            n <- length(x)
            structure(list(estimate=mean(x), sd = sqrt((n-1)/n)*sd(x)),
                      class = "fitdistr")
        }
        if(distname == "weibull" && is.null(start)) {
            ## log-Weibull is Gumbel, so start from that
            m <- mean(log(x)); v <- var(log(x))
            shape <- 1.2/sqrt(v); scale <- exp(m + 0.572/shape)
            start <- list(shape = shape, scale = scale)
            start <- start[!is.element(names(start), dots)]
        }
        if(distname == "gamma" && is.null(start)) {
            m <- mean(x); v <- var(x)
            start <- list(shape = m^2/v, rate = m/v)
            start <- start[!is.element(names(start), dots)]
        }
        if(distname == "uniform" && is.null(start)) {
            start <- list(min = min(x), max = max(x))
            start <- start[!is.element(names(start), dots)]
        }
        if(distname == "negative binomial" && is.null(start)) {
            m <- mean(x); v <- var(x)
            size <- if(v > m) m^2/(v - m) else 100
            start <- list(size = size, mu = m)
            start <- start[!is.element(names(start), dots)]
        }
        if(is.element(distname, c("cauchy", "logistic")) && is.null(start)) {
            start <- list(location = median(x), scale = IQR(x)/2)
            start <- start[!is.element(names(start), dots)]
        }
        if(distname == "t" && is.null(start)) {
            start <- list(m = median(x), s = IQR(x)/2, df = 10)
            start <- start[!is.element(names(start), dots)]
        }
    }
    if(is.null(start) || !is.list(start))
        stop("'start' must be a named list")
    nm <- names(start)
    ## reorder arguments to densfun
    f <- formals(densfun)
    args <- names(f)
    m <- match(nm, args)
    if(any(is.na(m)))
        stop("'start' specifies names which are not arguments to 'densfun'")
    formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
    dens <- function(parm, x, ...) densfun(x, parm, ...)
    if((l <- length(nm)) > 1)
        body(dens) <-
            parse(text = paste("densfun(x,",
                  paste("parm[", 1:l, "]", collapse = ", "),
                  ", ...)"))
    if("log" %in% args)
        res <- optim(start, mylogfn, x = x, hessian = TRUE, ...)
    else
        res <- optim(start, myfn, x = x, hessian = TRUE, ...)
    if(res$convergence > 0) stop("optimization failed")
    sds <- sqrt(diag(solve(res$hessian)))
    structure(list(estimate = res$par, sd = sds), class = "fitdistr")
}

print.fitdistr <-
    function(x, digits = getOption("digits"), ...)
{
    ans <- format(rbind(x$estimate, x$sd), digits=digits)
    ans[1, ] <- sapply(ans[1, ], function(x) paste("", x))
    ans[2, ] <- sapply(ans[2, ], function(x) paste("(", x, ")", sep=""))
    dn <- dimnames(ans)
    dn[[1]] <- rep("", 2)
    dn[[2]] <- paste(substring("      ", 1, (nchar(ans[2,]) - nchar(dn[[2]])) %/% 2), dn[[2]])
    dn[[2]] <- paste(dn[[2]], substring("      ", 1, (nchar(ans[2,]) - nchar(dn[[2]])) %/% 2))
    dimnames(ans) <- dn
    print(ans, quote = FALSE)
    x
}

coef.fitdistr <- function(object, ...) object$estimate

