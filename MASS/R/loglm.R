# file MASS/loglm.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
denumerate <- function(x) UseMethod("denumerate")

renumerate <- function(x) UseMethod("renumerate")

denumerate.formula <- function(x)
{
    if(length(x) == 1) {
        if(mode(x) == "numeric" ||
           (mode(x) == "name" &&
            any(substring(as.character(x), 1, 1) == as.character(1:9))))
            x <- as.name(paste(".v", x, sep = ""))
    } else {
        x[[2]] <- Recall(x[[2]])
        if(length(x) == 3 && x[[1]] != as.name("^"))
            x[[3]] <- Recall(x[[3]])
    }
    x
}

renumerate.formula <- function(x)
{
    if(length(x) == 1) {
        if(mode(x) == "name"
           && nchar(xx <- as.character(x)) > 2
           && substring(xx, 1, 2) == ".v")
            x <- as.name(substring(xx, 3))
    } else {
        x[[2]] <- Recall(x[[2]])
        if(length(x) == 3 && x[[1]] != as.name("^"))
           x[[3]] <- Recall(x[[3]])
    }
    x
}

loglm <-
  function(formula, data, subset, na.action, ...)
{
    assign(".call", match.call(), envir=.GlobalEnv)
    if(missing(data) || inherits(data, "data.frame")) {
        m <- match.call(expand = FALSE)
        m$... <- NULL
        m[[1]] <- as.name("model.frame")
        data <- eval.parent(m)
        assign(".formula", as.formula(attr(data, "terms")), envir=.GlobalEnv)
    } else {
        trms <- attr(data, "terms") <- terms(formula <- denumerate(formula))
        assign(".formula", renumerate(as.formula(trms)), envir=.GlobalEnv)
    }
    loglm1(formula, data, ...)
}

loglm1 <- function(formula, data, ...) UseMethod("loglm1", data)

loglm1.xtabs <-
function(formula, data, ...)
{
    attr(data, "marginals") <- attr(data, "call") <- class(data) <- NULL
    NextMethod("loglm1")
}

loglm1.data.frame <-
function(formula, data, ...)
{
    trms <- attr(data, "terms")
    if(is.null(trms)) stop("'data' has no 'terms' attribute")
    if(attr(trms, "response") == 0) stop("formula specifies no response")
    resp <- match(as.character(attr(trms, "variables"))[1+attr(trms, "response")],
                  names(data))
    off <- attr(trms, "offset")
    fac <- data.frame(lapply(data[, -c(resp, off)], as.factor))
    rsp <- data[, resp]
    tab <- table(fac)
    if(max(tab) > 1) {
#
# an extra factor needed for repeated frequencies
#
        i <- do.call("order", rev(fac))
        fac <- fac[i,  ]
        rsp <- rsp[i]
        fac$.Within. <-
            factor(unlist(sapply(tab,
                                 function(x) if(x > 0) seq(x) else NULL)))
    }
    dn <- lapply(fac, levels)
    dm <- sapply(dn, length)
    offset <- model.offset(data)
    if (is.null(offset)) offset <- 0
    offset <- rep(offset, length.out = nrow(data))
    data <- structure(array(-1, dm, dn), terms = trms)
    data[do.call("cbind", lapply(fac, as.numeric))] <- rsp
    st <- array(as.numeric(data >= 0), dm, dn)
    st[do.call("cbind", lapply(fac, as.numeric))] <- exp(offset)
    data[data < 0] <- 0
    loglm1.default(formula, data, ..., start = st)
}

loglm1.default <-
function(formula, data, start = rep(1, length(data)), fitted = FALSE,
	keep.frequencies = fitted, param = TRUE, eps =
	1/10, iter = 40, print = FALSE, ...)
{
    trms <- attr(data, "terms")
    if(is.null(trms)) stop("'data' has no 'terms' attribute")
    factors <- attr(trms, "factors") > 0
    if((r <- attr(trms, "response")))
        factors <- factors[-r,  , drop = FALSE]
    nt <- ncol(factors)
    fo <- order(colSums(factors))
    factors <- factors[, fo, drop = FALSE]
    ff <- crossprod(factors)
    keep <- rep(TRUE, nt)
    j <- 0
    while((j <- j + 1) < nt) keep[j] <- ff[j, j] > max(ff[j, (j + 1):nt])
    factors <- factors[, keep, drop = FALSE]
    ldim <- length(dim(data))
    nnames <- paste(".v", 1:ldim, sep = "")
    which <- structure(1:ldim, names = nnames)
    if(!is.null(anames <- names(dimnames(data))))
        which <- c(which, structure(which, names = anames))
    margins <- apply(factors, 2, function(x, which, nam)
                     as.vector(which[nam[x]]), which, rownames(factors))
    if(is.matrix(margins))
        margins <- as.list(data.frame(margins))
    else margins <- structure(as.list(margins), names = names(margins))
    Fit <- loglin(data, margins, start = start, fit = fitted,
                  param = param, eps = eps, iter = iter, print = print)
    if(exists(".formula")) {
        Fit$call <- .call
        Fit$formula <- .formula
    }
    class(Fit) <- "loglm"
    if(keep.frequencies) Fit$frequencies <- structure(data, terms = NULL)
    if(fitted) {
        names(Fit)[match("fit", names(Fit))] <- "fitted"
        attr(Fit$fitted, "terms") <- NULL
    }
    Fit$deviance <- Fit$lrt
    Fit$nobs <- length(data)
    Fit$df <- Fit$df - sum(start == 0)
    Fit$terms <- trms # for stepAIC
    Fit
}


anova.loglm <- function(object, ..., test = c("Chisq", "chisq", "LR"))
{
    test <- match.arg(test)
    margs <- function(...) nargs()
    if(!(k <- margs(...))) return(object)
    objs <- list(object, ...)
    dfs <- sapply(objs, "[[", "df")
    o <- order( - dfs)
    objs <- objs[o]
    dfs <- c(dfs[o], 0)
    forms <- lapply(objs, formula)
    dev <- c(sapply(objs, "[[", "lrt"), 0)
    M <- array(0, c(k + 2, 5),
               list(c(paste("Model", 1:(k + 1)), "Saturated"),
                    c("Deviance", "df", "Delta(Dev)", "Delta(df)", "P(> Delta(Dev)")))
    M[, 1] <- dev
    M[, 2] <- dfs
    M[-1, 3] <- dev[1:(k + 1)] - dev[2:(k + 2)]
    M[-1, 4] <- dfs[1:(k + 1)] - dfs[2:(k + 2)]
    M[-1, 5] <- 1 - pchisq(M[-1, 3], M[-1, 4])
    res <- structure(M, formulae = forms)
    class(res) <- "anova.loglm"
    res
}

print.anova.loglm <- function(x, ...)
{
    rjustify <- function(str) {
        m <- max(n <- nchar(str))
        blanks <- format(c("", str[n == m][1]))[1]
        paste(substring(blanks, 0, m - n), str, sep = "")
    }
    y <- x
    y[, 5] <- round(y[, 5], 5)
    R <- array("", dim(x), dimnames(x))
    for(j in 1:5) {
        colj <- rjustify(c(colnames(x)[j], format(y[, j])))
        R[, j] <- colj[-1]
        colnames(R)[j] <- colj[1]
    }
    R[1, 3:5] <- ""
    pform <- function(form)
        if(length(form) == 2) form else form[c(2, 1, 3)]
    forms <- attr(x, "formulae")
    cat("LR tests for hierarchical log-linear models\n\n")
    for(i in seq_along(forms))
        cat(paste("Model ", i, ":\n", sep = ""),
            deparse(pform(forms[[i]])), "\n")
    cat("\n")
    print(R, quote = FALSE)
    invisible(x)
}

print.loglm <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    ts.array <- rbind(c(x$lrt, x$df,
                        if(x$df > 0) 1 - pchisq(x$lrt, x$df) else 1),
                      c(x$pearson, x$df,
                        if(x$df > 0) 1 - pchisq(x$pearson, x$df)
                        else 1))
    dimnames(ts.array) <- list(c("Likelihood Ratio",
                                 "Pearson"), c("X^2", "df", "P(> X^2)"))
    cat("\nStatistics:\n")
    print(ts.array)
    invisible(x)
}

summary.loglm <- function(object, fitted = FALSE, ...)
{
    ts.array <- rbind(c(object$lrt, object$df,
                        if(object$df > 0) 1 - pchisq(object$lrt, object$df)
                        else 1), c(object$pearson, object$df,
                                   if(object$df > 0)
                                   1 - pchisq(object$pearson, object$df)
                                   else 1))
    dimnames(ts.array) <- list(c("Likelihood Ratio", "Pearson"),
                               c("X^2", "df", "P(> X^2)"))
    if(fitted) {
        if(is.null(object$fitted) || is.null(object$freqencies)) {
            cat("Re-fitting to find fitted values\n")
            object <- update(object, fitted = TRUE, keep.frequencies = TRUE)
        }
        fit <- format(round(object$fit, 1))
        OE <- array(paste(format(object$freq), " (", fit, ")", sep = ""),
                    dim(fit), dimnames(object$freq))
    }  else OE <- NULL
    structure(list(formula = formula(object), tests = ts.array, oe = OE),
              class = "summary.loglm")
}

print.summary.loglm <- function(x, ...)
{
    cat("Formula:\n")
    print(formula(x))
    cat("\nStatistics:\n")
    print(x$tests)
    if(!is.null(x$oe)) {
        cat("\nObserved (Expected):\n")
        print(x$oe, quote = FALSE)
    }
    invisible(x)
}

update.loglm <- function (object, formula, ...)
{
    if (is.null(call <- object$call))
        stop("'object' has no 'call' component.  Updating not possible")
    if (fix <- !missing(formula)) {
        object$formula <- denumerate(object$formula)
        formula <- denumerate(as.formula(formula))
        call$formula <- update.formula(formula(object), formula)
    }
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        ## do these individually to allow NULL to remove entries.
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    result <- eval.parent(call)
    if (fix) {
        form <- renumerate(result$formula)
        result$call$formula <- unclass(result$formula <- form)
    }
    result
}

fitted.loglm <- function(object, ...)
{
    if(!is.null(object$fit))
        return(unclass(object$fit))
    cat("Re-fitting to get fitted values\n")
    unclass(update(object, fitted = TRUE, keep.frequencies = FALSE)$fitted)
}

residuals.loglm <-
    function(object, type = c("deviance", "pearson", "response"), ...)
{
    type <- match.arg(type)
    if(is.null(object$fit) || is.null(object$freq)) {
        cat("Re-fitting to get frequencies and fitted values\n")
        object <- update(object, fitted = TRUE, keep.frequencies = TRUE)
    }
    y <- object$freq
    mu <- object$fit
    res <- y - mu
    nz <- mu > 0
    y <- y[nz]
    mu <- mu[nz]
    res[nz] <-
        switch(type,
               deviance = sign(y - mu) *
                 sqrt(2*abs(y*log((y + (y == 0))/mu) - (y - mu))),
               pearson = (y - mu)/sqrt(mu),
               response = y - mu)
    res
}

coef.loglm <- function(object, ...)
{
    if(!is.null(cf <- object$param)) return(cf)
    cat("Re-fitting to calculate missing coefficients\n")
    update(object, param = TRUE)$param
}
