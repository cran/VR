glmmPQL <- function(fixed, random, family, data, correlation, weights,
                    control, niter = 10, verbose = TRUE, ...)
{
    if(!require("nlme")) stop("package 'nlme' is essential")
    ## family
    if(is.character(family)) family <- get(family)
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {
	print(family)
	stop("'family' not recognized")
    }
    m <- mcall <- Call <- match.call()
    nm <- names(m)[-1]
    keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
    for(i in nm[!keep]) m[[i]] <- NULL
    allvars <-
        if (is.list(random))
            allvars <- c(all.vars(fixed), names(random),
                         unlist(lapply(random, function(x) all.vars(formula(x)))))
        else c(all.vars(fixed), all.vars(random))
    ## allvars does not contain offset term.
    Terms <- terms(fixed, data=data)
    off <- attr(Terms, "offset")
    if(length(off<- attr(Terms, "offset"))) allvars <-
        c(allvars, as.character(attr(Terms, "variables"))[off+1])
    m$formula <- as.formula(paste("~", paste(allvars, collapse="+")))
    environment(m$formula) <- environment(fixed)
    m$drop.unused.levels <- TRUE
    m[[1]] <- as.name("model.frame")
    mf <- eval.parent(m)
    off <- model.offset(mf)
    if(is.null(off)) off <- 0
    w <-  model.weights(mf)
    if(is.null(w)) w <- rep(1, nrow(mf))
    mf$wts <- w
    fit0 <- glm(formula=fixed, family=family, data=mf, weights = wts, ...)
    w <- fit0$prior.weights
    eta <- fit0$linear.predictor
    zz <- eta + fit0$residuals - off
    wz <- fit0$weights
    fam <- family

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
        eta <- fitted(fit) + off
        if(sum((eta-etaold)^2) < 1e-6*sum(eta^2)) break;
        mu <- fam$linkinv(eta)
        mu.eta.val <- fam$mu.eta(eta)
        mf$zz <- eta + (fit0$y - mu)/mu.eta.val - off
        wz <- w * mu.eta.val^2 / fam$variance(mu)
        mf$invwt <- 1/wz
        mcall$data <- mf
    }
    attributes(fit$logLik) <- NULL # needed for some versions of nlme
    fit$call <- Call
    fit$family <- family
    oldClass(fit) <- c("glmmPQL", oldClass(fit))
    fit
}

predict.glmmPQL <-
  function(object, newdata = NULL, type = c("link", "response"),
           level = Q, na.action = na.pass, ...)
{
    type <- match.arg(type)
    Q <- object$dims$Q
    if(missing(newdata)) {
        pred <- fitted(object, level = level)
        pred <- switch(type,
                       link = pred,
                       response = object$family$linkinv(pred))
        if(!is.null(object$na.action))
            pred <- napredict(object$na.action, pred)
    } else {
        class(object) <- class(object)[-1]
        pred <- predict(object, newdata, level = level,
                        na.action = na.action)
        switch(type,
               response = {pred <- object$family$linkinv(pred)},
               link =)
    }
    pred
}
