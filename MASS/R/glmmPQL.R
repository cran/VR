glmmPQL <- function(fixed, random, family, data, correlation, weights,
                    control,
                    niter = 10, verbose = TRUE, ...)
{
    ## family
    if(is.character(family)) family <- get(family)
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {
	print(family)
	stop("`family' not recognized")
    }
    m <- mcall <- Call <- match.call()
    nm <- names(m)[-1]
    keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
    for(i in nm[!keep]) m[[i]] <- NULL
    allvars <- c(all.vars(fixed), all.vars(random))
    m$formula <- as.formula(paste("~", paste(allvars, collapse="+")))
    m$drop.unused.levels <- T
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    off <- model.extract(mf, "offset")
    if(is.null(off)) off <- 0
    w <-  model.extract(mf, "weights")
    if(is.null(w)) w <- rep(1, nrow(mf))
    mf$wts <- w
    fit0 <- glm(formula=fixed, family=family, data=mf, weights = wts, ...)
    w <- fit0$prior.weights
    eta <- fit0$linear.predictor - off
    zz <- eta + fit0$residuals
    wz <- fit0$weights
    fam <- family
    ## na.action fix here
    nm <- names(mcall)[-1]
    keep <- is.element(nm, c("fixed", "random", "data", "subset",
                             "na.action", "control"))
    for(i in nm[!keep]) mcall[[i]] <- NULL
    fixed[[2]] <- quote(zz)
    mcall[["fixed"]] <- fixed
    mcall[[1]] <- as.name("lme")
    mcall$random <- random
    mcall$method <- "ML"
    if(!missing(correlation))
        mcall$correlation <- correlation
    mcall$weights <- quote(varFixed(~invwt))
    mf$zz <- zz
    mf$invwt <- 1/wz
    mcall$data <- mf
    for(i in 1:niter) {
        if(verbose) cat("iteration", i, "\n")
        fit <- eval(mcall)
        etaold <- eta
        ##update zz and invwt
        eta <- fit$fitted[, 2]
        if(sum((eta-etaold)^2) < 1e-6*sum(eta^2)) break;
        mu <- fam$linkinv(eta)
        mu.eta.val <- fam$mu.eta(eta)
        mf$zz <- eta + (fit0$y - mu)/mu.eta.val  - off
        wz <- w * mu.eta.val^2 / fam$variance(mu)
        mf$invwt <- 1/wz
        mcall$data <- mf
    }
    attributes(fit$logLik) <- NULL # needed for some versions of nlme
    fit$call <- Call
    fit
}

