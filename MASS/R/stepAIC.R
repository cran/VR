# file MASS/stepAIC.q
# copyright (C) 1994-2003 W. N. Venables and B. D. Ripley
#
stepAIC <-
  function(object, scope, scale = 0,
           direction = c("both", "backward", "forward"),
           trace = 1, keep = NULL, steps = 1000, use.start = FALSE, k = 2, ...)
{
    mydeviance <- function(x, ...)
    {
        dev <- deviance(x)
        if(!is.null(dev)) dev else extractAIC(x, k=0)[2]
    }

    cut.string <- function(string)
    {
        if(length(string) > 1)
            string[-1] <- paste("\n", string[-1], sep = "")
        string
    }

    re.arrange <- function(keep)
    {
        namr <- names(k1 <- keep[[1]])
        namc <- names(keep)
        nc <- length(keep)
        nr <- length(k1)
        array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, namc))
    }

    step.results <- function(models, fit, object, usingCp=FALSE)
    {
        change <- sapply(models, "[[", "change")
        rd <- sapply(models, "[[", "deviance")
        dd <- c(NA, abs(diff(rd)))
        rdf <- sapply(models, "[[", "df.resid")
        ddf <- c(NA, abs(diff(rdf)))
        AIC <- sapply(models, "[[", "AIC")
        heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                     "\nInitial Model:", deparse(as.vector(formula(object))),
                     "\nFinal Model:", deparse(as.vector(formula(fit))),
                     "\n")
        aod <-
            if(usingCp)
                data.frame(Step = change, Df = ddf, Deviance = dd,
                           "Resid. Df" = rdf, "Resid. Dev" = rd,
                           Cp = AIC, check.names = FALSE)
            else data.frame(Step = change, Df = ddf, Deviance = dd,
                            "Resid. Df" = rdf, "Resid. Dev" = rd,
                            AIC = AIC, check.names = FALSE)
        attr(aod, "heading") <- heading
        class(aod) <- c("Anova", "data.frame")
        fit$anova <- aod
        fit
    }

    Terms <- terms(object)
    object$formula <- Terms
    if(inherits(object, "lme")) object$call$fixed <- Terms
    else if(inherits(object, "gls")) object$call$model <- Terms
    else object$call$formula <- Terms
    if(use.start) warning("'use.start' cannot be used with R's version of glm")
    md <- missing(direction)
    direction <- match.arg(direction)
    backward <- direction == "both" | direction == "backward"
    forward <- direction == "both" | direction == "forward"
    if(missing(scope)) {
        fdrop <- numeric(0)
        fadd <- attr(Terms, "factors")
        if(md) forward <- FALSE
    } else {
        if(is.list(scope)) {
            fdrop <- if(!is.null(fdrop <- scope$lower))
                attr(terms(update.formula(object, fdrop)), "factors")
            else numeric(0)
            fadd <- if(!is.null(fadd <- scope$upper))
                attr(terms(update.formula(object, fadd)), "factors")
        } else {
            fadd <- if(!is.null(fadd <- scope))
                attr(terms(update.formula(object, scope)), "factors")
            fdrop <- numeric(0)
        }
    }
    models <- vector("list", steps)
    if(!is.null(keep)) keep.list <- vector("list", steps)
    ## watch out for partial matching here.
    if(is.list(object) && (nmm <- match("nobs", names(object), 0)) > 0)
        n <- object[[nmm]]
    else n <- length(residuals(object))
    fit <- object
    bAIC <- extractAIC(fit, scale, k = k, ...)
    edf <- bAIC[1]
    bAIC <- bAIC[2]
    if(is.na(bAIC))
        stop("AIC is not defined for this model, so stepAIC cannot proceed")
    nm <- 1
    Terms <- terms(fit)
    if(trace) {
        cat("Start:  AIC=", format(round(bAIC, 2)), "\n",
            cut.string(deparse(as.vector(formula(fit)))), "\n\n")
	flush.console()
    }
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - edf,
                         change = "", AIC = bAIC)
    if(!is.null(keep)) keep.list[[nm]] <- keep(fit, bAIC)
    usingCp <- FALSE
    while(steps > 0) {
        steps <- steps - 1
        AIC <- bAIC
        ffac <- attr(Terms, "factors")
        ## don't drop strata terms
        if(!is.null(sp <- attr(Terms, "specials")) &&
           !is.null(st <- sp$strata)) ffac <- ffac[-st,]
        scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
        aod <- NULL
        change <- NULL
        if(backward && length(scope$drop)) {
            aod <- dropterm(fit, scope$drop, scale = scale,
                            trace = max(0, trace - 1), k = k, ...)
            rn <- row.names(aod)
            row.names(aod) <- c(rn[1], paste("-", rn[-1], sep=" "))
            ## drop all zero df terms first.
            if(any(aod$Df == 0, na.rm=TRUE)) {
                zdf <- aod$Df == 0 & !is.na(aod$Df)
                nc <- match(c("Cp", "AIC"), names(aod))
                nc <- nc[!is.na(nc)][1]
                ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
                if(any(ch)) {
                    warning("0 df terms are changing AIC")
                    zdf <- zdf[!ch]
                }
                ## drop zero df terms first: one at time since they
                ## may mask each other
                if(length(zdf) > 0)
                    change <- rev(rownames(aod)[zdf])[1]
            }
        }
        if(is.null(change)) {
            if(forward && length(scope$add)) {
                aodf <- addterm(fit, scope$add, scale = scale,
                                trace = max(0, trace - 1), k = k, ...)
                rn <- row.names(aodf)
                row.names(aodf) <- c(rn[1], paste("+", rn[-1], sep=" "))
                aod <-
                    if(is.null(aod)) aodf
                    else rbind(aod, aodf[-1, , drop=FALSE])
            }
            attr(aod, "heading") <- NULL
            if(is.null(aod) || ncol(aod) == 0) break
            ## need to remove any terms with zero df from consideration
            nzdf <- if(!is.null(aod$Df)) aod$Df != 0 | is.na(aod$Df)
            aod <- aod[nzdf, ]
            if(is.null(aod) || ncol(aod) == 0) break
            nc <- match(c("Cp", "AIC"), names(aod))
            nc <- nc[!is.na(nc)][1]
            o <- order(aod[, nc])
            if(trace) {
		print(aod[o,  ])
		flush.console()
	    }
            if(o[1] == 1) break
            change <- rownames(aod)[o[1]]
        }
        usingCp <- match("Cp", names(aod), 0) > 0
        ## may need to look for a 'data' argument in parent
	fit <- update(fit, paste("~ .", change), evaluate = FALSE)
        fit <- eval.parent(fit)
        if(is.list(fit) && (nmm <- match("nobs", names(fit), 0)) > 0)
            nnew <- fit[[nmm]]
        else nnew <- length(residuals(fit))
        if(nnew != n)
            stop("number of rows in use has changed: remove missing values?")
        Terms <- terms(fit)
        bAIC <- extractAIC(fit, scale, k = k, ...)
        edf <- bAIC[1]
        bAIC <- bAIC[2]
        if(trace) {
            cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n",
                cut.string(deparse(as.vector(formula(fit)))), "\n\n")
	    flush.console()
	}
        ## add a tolerance as dropping 0-df terms might increase AIC slightly
        if(bAIC >= AIC + 1e-7) break
        nm <- nm + 1
        models[[nm]] <-
            list(deviance = mydeviance(fit), df.resid = n - edf,
                 change = change, AIC = bAIC)
        if(!is.null(keep)) keep.list[[nm]] <- keep(fit, bAIC)
    }
    if(!is.null(keep)) fit$keep <- re.arrange(keep.list[seq(nm)])
    step.results(models = models[seq(nm)], fit, object, usingCp)
}

extractAIC.loglm <- function(fit, scale, k = 2, ...)
{
    edf <- fit$n - fit$df
    c(edf,  fit$deviance + k * edf)
}

extractAIC.lme <- function(fit, scale, k = 2, ...)
{
    if(fit$method != "ML") stop("AIC undefined for REML fit")
    res <- logLik(fit)
    edf <- attr(res, "df")
    c(edf,  -2*res + k * edf)
}

extractAIC.gls <- function(fit, scale, k = 2, ...)
{
    if(fit$method != "ML") stop("AIC undefined for REML fit")
    res <- logLik(fit)
    edf <- attr(res, "df")
    c(edf,  -2*res + k * edf)
}

terms.gls <- terms.lme <- function(x, ...) terms(formula(x), ...)
