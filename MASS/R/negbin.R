# file MASS/negbin.q
# copyright (C) 1994-2000 W. N. Venables and B. D. Ripley
#
"anova.negbin"<-
function(object, ..., test = "Chisq")
{
  dots <- list(...)
  if(length(dots) == 0) {
    warning("tests made without re-estimating theta")
    object$call[[1]] <- as.name("glm")
    if(is.null(object$link))
      object$link <- as.name("log")
    object$call$family <- call("negative.binomial", theta = object$
                               theta, link = object$link)
    anova.glm(object, test = test)
  } else {
    if(test != "Chisq")
      warning("only Chi-squared LR tests are implemented")
    mlist <- list(object, ...)
    nt <- length(mlist)
    dflis <- sapply(mlist, function(x) x$df.resid)
    s <- sort.list(-dflis)
    mlist <- mlist[s]
    if(any(!sapply(mlist, inherits, "negbin")))
      stop("not all objects are of class `negbin'")
    rsp <- unique(sapply(mlist, function(x) paste(formula(x)[2])))
    mds <- sapply(mlist, function(x) paste(formula(x)[3]))
    ths <- sapply(mlist, function(x) x$theta)
    dfs <- dflis[s]
    lls <- sapply(mlist, function(x) x$twologlik)
    tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
    df <- c(NA,  - diff(dfs))
    x2 <- c(NA, diff(lls))
    pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
    out <- data.frame(Model = mds, theta = ths, Resid.df = dfs,
                      "2 x log-lik." = lls, Test = tss, df = df, LRtest = x2,
                      Prob = pr)
    names(out) <- c("Model", "theta", "Resid. df",
                    "   2 x log-lik.", "Test", "   df", "LR stat.", "Pr(Chi)")
    class(out) <- c("Anova", "data.frame")
    attr(out, "heading") <-
      c("Likelihood ratio tests of Negative Binomial Models\n",
        paste("Response:", rsp))
    out
  }
}

print.Anova <- function(x, ...)
{
    heading <- attr(x, "heading")
    if(!is.null(heading)) cat(heading, sep = "\n")
    attr(x, "heading") <- NULL
    print.data.frame(x)
}

"family.negbin"<-
function(object, ...)
    object$family

"glm.convert"<-
function(object)
{
    object$call[[1]] <- as.name("glm")
    if(is.null(object$link))
        object$link <- as.name("log")
    object$call$family <- call("negative.binomial", theta = object$theta,
                               link = object$link)
    object$call$init.theta <- object$call$link <- NULL
    class(object) <- c("glm", "lm")
    object
}

glm.nb <- function(formula = formula(data), data = parent.frame(), weights,
		   subset, na.action, start = NULL, etastart = NULL,
		   control = glm.control(...), method = "glm.fit",
		   model = TRUE,
                   x = FALSE, y = TRUE, contrasts = NULL, ...,
		   init.theta, link = log)
{
    loglik <- function(n, th, mu, y)
    {
        sum(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
            y * log(mu + (y == 0)) - (th + y) * log(th + mu))
    }
    link <- substitute(link)
    if(missing(init.theta)) {
        fam0 <- do.call("poisson", list(link = link))
    } else {
        fam0 <- do.call("negative.binomial", list(theta = init.theta, link = link))
    }
    Call <- match.call()
    m <- match.call(expand = FALSE)
    m$method <- m$model <- m$x <- m$y <- m$control <- m$contrasts <-
        m$init.theta <- m$link <- m$start <- m$etastart <- m$... <-  NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    if(method == "model.frame") return(m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if(length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev <- xlev[!sapply(xlev, is.null)]
    } else xlev <- NULL
    a <- attributes(m)
    Y <- model.extract(m, response)
    X <- model.matrix(Terms, m, contrasts)
    w <- model.extract(m, weights)
    if(!length(w)) w <- rep(1, nrow(m))
    else if(any(w < 0)) stop("negative weights not allowed")
    offset <- model.extract(m, offset)
    n <- length(Y)
    if(!is.null(method)) {
        if(!exists(method, mode = "function"))
            stop(paste("unimplemented method:", method))
    }
    else method <- "glm.fit"
    glm.fitter <- get(method)
    if(control$trace > 1) cat("Initial fit:\n")
    fit <- glm.fitter(x = X, y = Y, w = w, start = start,
                      etastart = etastart,
                      offset = offset, family = fam0,
                      control = list(maxit=control$maxit,
                      epsilon = control$epsilon,
                      trace = control$trace > 1))
    class(fit) <- c("glm", "lm")
    mu <- fit$fitted
    th <- as.vector(theta.ml(Y, mu, n, limit=control$maxit, trace =
                             control$trace> 2))
    if(control$trace > 1)
        cat("Initial value for theta:", signif(th), "\n")
    fam <- do.call("negative.binomial", list(theta = th, link = link))
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    g <- fam$linkfun
    Lm <- loglik(n, th, mu, Y)
    Lm0 <- Lm + 2 * d1
    while((iter <- iter + 1) <= control$maxit &&
          (abs(Lm0 - Lm)/d1 + abs(del)/d2) > control$epsilon) {
        eta <- g(mu)
        fit <- glm.fitter(x = X, y = Y, w = w, etastart =
                          eta, offset = offset, family = fam,
                          control = list(maxit=control$maxit,
                          epsilon = control$epsilon,
                          trace = control$trace > 1),
                          intercept = attr(Terms, "intercept") > 0)
        t0 <- th
        th <- theta.ml(Y, mu, n, limit=control$maxit, trace = control$trace > 2)
        fam <- do.call("negative.binomial", list(theta = th, link = link))
        mu <- fit$fitted
        del <- t0 - th
        Lm0 <- Lm
        Lm <- loglik(n, th, mu, Y)
        if(control$trace) {
            Ls <- loglik(n, th, Y, Y)
            Dev <- 2 * (Ls - Lm)
            cat("Theta(", iter, ") =", signif(th),
                ", 2(Ls - Lm) =", signif(Dev), "\n")
        }
    }
    if(!is.null(attr(th, "warn"))) fit$th.warn <- attr(th, "warn")
    if(iter > control$maxit) {
        warning("alternation limit reached")
        fit$th.warn <- "alternation limit reached"
    }

  # If an offset and intercept are present, iterations are needed to
  # compute the Null deviance; these are done here, unless the model
  # is NULL, in which case the computations have been done already
  #
    if(any(offset) && attr(Terms, "intercept")) {
        null.deviance <-
            if(length(Terms))
                glm.fitter(X[, "(Intercept)", drop = FALSE], Y, w,
                           offset = offset, family = fam,
                           control = list(maxit=control$maxit,
                           epsilon = control$epsilon, trace = control$trace > 1)
                           )$deviance
           else fit$deviance
        fit$null.deviance <- null.deviance
    }
    class(fit) <- c("negbin", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    Call$init.theta <- as.vector(th)
    Call$link <- link
    fit$call <- Call
    if(model) fit$model <- m
    if(x) fit$x <- X
    if(!y) fit$y <- NULL
    fit$theta <- as.vector(th)
    fit$SE.theta <- attr(th, "SE")
    fit$twologlik <- as.vector(2 * Lm)
    fit$aic <- -fit$twologlik + 2*fit$rank + 2
    fit$contrasts <- if (length(clv <- unlist(lapply(m, class))))
        options("contrasts")[[1]] else FALSE
    fit$xlevels <- xlev
    fit$method <- method
    fit$control <- control
    fit
}

"negative.binomial" <-
function(theta = stop("theta must be specified"), link = "log")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp)) {
        linktemp <- deparse(linktemp)
        if (linktemp == "link")
            linktemp <- eval(link)
    }
    if (any(linktemp == c("log", "identity", "sqrt")))
        stats <- make.link(linktemp)
    else stop(paste(linktemp, "link not available for negative binomial",
                    "family; available links are", "\"identity\", \"log\" and \"sqrt\""))
    .Theta <- theta
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
    famname <- paste("Negative Binomial(", format(round(theta, 4)), ")",
                     sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                   linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                   aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                   validmu = validmu, valideta = stats$valideta), class = "family")
}

rnegbin <-
function(n, mu = n, theta = stop("theta must be given"))
{
    k <- if(length(n) > 1) length(n) else n
    rpois(k, (mu * rgamma(k, theta))/theta)
}

"summary.negbin" <-
function(object, dispersion = 1, correlation = TRUE, ...)
{
    if(is.null(dispersion)) dispersion <- 1
    summ <- c(summary.glm(object, dispersion = dispersion,
                          correlation = correlation),
              object[c("theta", "SE.theta", "twologlik", "th.warn")])
    class(summ) <- c("summary.negbin", "summary.glm")
    summ
}

"print.summary.negbin" <- function(x, ...)
{
    NextMethod()
    dp <- 2 - floor(log10(x$SE.theta))
    cat("\n              Theta: ", format(round(x$theta, dp), nsmall=dp),
        "\n          Std. Err.: ", format(round(x$SE.theta, dp), nsmall=dp),
        "\n")
    if(!is.null(x$th.warn))
        cat("Warning while fitting theta:", x$th.warn,"\n")
    cat("\n 2 x log-likelihood: ", format(round(x$twologlik, 3), nsmall=dp), "\n")
    invisible(x)
}

"theta.md" <-
    function(y, u, dfr, limit = 20, eps = .Machine$double.eps^0.25)
{
    if(inherits(y, "lm")) {
        u <- y$fitted
        dfr <- y$df.residual
        y <- if(is.null(y$y)) u + residuals(y) else y$y
    }
    n <- length(y)
    t0 <- n/sum((y/u - 1)^2)
    a <- 2 * sum(y * log(pmax(1, y)/u)) - dfr
    it <- 0
    del <- 1
    while((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        top <- a - 2 * sum((y + t0) * log((y + t0)/(u + t0)))
        bot <- 2 * sum((y - u)/(u + t0) - log((y + t0)/(u + t0)))
        del <- top/bot
        t0 <- t0 - del
    }
    if(t0 < 0) {
        t0 <- 0
        warning("estimator truncated at zero")
        attr(t0, "warn") <- "estimate truncated at zero"
    }
    t0
}

"theta.ml" <-
    function(y, mu, n = length(y), limit = 10, eps = .Machine$double.eps^0.25,
             trace=FALSE)
{
    score <- function(n, th, mu, y)
        sum(digamma(th + y) - digamma(th) + log(th) +
            1 - log(th + mu) - (y + th)/(mu + th))
    info <- function(n, th, mu, y)
        sum( - trigamma(th + y) + trigamma(th) - 1/th +
            2/(mu + th) - (y + th)/(mu + th)^2)
    if(inherits(y, "lm")) {
        mu <- y$fitted
        y <- if(is.null(y$y)) mu + residuals(y) else y$y
    }
    t0 <- n/sum((y/mu - 1)^2)
    it <- 0
    del <- 1
    if(trace) cat("theta.ml: initial theta =", signif(t0), "\n")
    while((it <- it + 1) < limit && abs(del) > eps) {
        t0 <- abs(t0)
        del <- score(n, t0, mu, y)/(i <- info(n, t0, mu, y))
        t0 <- t0 + del
        if(trace) cat("theta.ml: iter", it," theta =", signif(t0), "\n")
    }
    if(t0 < 0) {
        t0 <- 0
        warning("estimator truncated at zero")
        attr(t0, "warn") <- "estimate truncated at zero"
    }
    if(it == limit) {
        warning("iteration limit reached")
        attr(t0, "warn") <- "iteration limit reached"
    }
    attr(t0, "SE") <- sqrt(1/i)
    t0
}

"theta.mm" <-
function(y, u, dfr, limit = 10, eps = .Machine$double.eps^0.25)
{
  if(inherits(y, "lm")) {
    u <- y$fitted
    dfr <- y$df.residual
    y <- if(is.null(y$y)) u + residuals(y) else y$y
  }
  n <- length(y)
  t0 <- n/sum((y/u - 1)^2)
  it <- 0
  del <- 1
  while((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    del <- (sum((y - u)^2/(u + u^2/t0)) - dfr)/sum((y - u)^2/(u + t0)^2)
    t0 <- t0 - del
  }
  if(t0 < 0) {
    t0 <- 0
    warning("estimator truncated at zero")
    attr(t0, "warn") <- "estimate truncated at zero"
  }
  t0
}

logLik.negbin <- function(object, ...)
{
    if (length(list(...)))
        warning("extra arguments discarded")
    p <- object$rank + 1 # for theta
    val <- object$twologlik/2
    attr(val, "df") <- p
    class(val) <- "logLik"
    val

}
