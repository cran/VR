# file MASS/lda.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
lda <- function(x, ...)
{
    if(is.null(class(x))) class(x) <- data.class(x)
    UseMethod("lda", x, ...)
}

lda.formula <- function(formula, data = NULL, ...,
			subset, na.action = na.fail)
{
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, sys.frame(sys.parent()))))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))
    Terms <- attr(m, "terms")
    grouping <- model.extract(m, "response")
    x <- model.matrix(Terms, m)
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0) x <- x[, -xint, drop=FALSE]
    res <- lda.default(x, grouping, ...)
    res$terms <- Terms
    res$call <- match.call()
    res
}

lda.data.frame <- function(x, ...)
{
    res <- lda.matrix(data.matrix(x), ...)
    res$call <- match.call()
    res
}

lda.Matrix <- function(x, ...)
{
    res <- lda.matrix(as.matrix(x), ...)
    res$call <- match.call()
    res
}

lda.matrix <- function(x, grouping, ...,
		       subset, na.action = na.fail)
{
    if(!missing(subset)) {
        x <- x[subset, , drop = FALSE]
        grouping <- grouping[subset]
    }
    if(!missing(na.action)) {
        dfr <- na.action(structure(list(g = grouping, x = x),
                                   class = "data.frame"))
        grouping <- dfr$g
        x <- dfr$x
    }
    res <- NextMethod("lda")
    res$call <- match.call()
    res
}

lda.default <-
  function(x, grouping, prior = proportions, tol = 1.0e-4,
           method = c("moment", "mle", "mve", "t"), CV=FALSE,
           nu = 5, ...)
{
    which.is.max <- function(x) {
        d <- (1:length(x))[x == max(x)]
        if(length(d) > 1) d <- sample(d, 1)
        d
    }

    if(is.null(dim(x))) stop("x is not a matrix")
    n <- nrow(x)
    p <- ncol(x)
    if(n != length(grouping)) stop("nrow(x) and length(grouping) are different")
    g <- as.factor(grouping)
    lev <- lev1 <- levels(g)
    counts <- as.vector(table(g))
    if(any(counts == 0)) {
        warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),
                      "are empty"))
        lev1 <- lev[counts > 0]
        g <- factor(g, levels=lev1)
        counts <- as.vector(table(g))
    }
    proportions <- counts/n
    ng <- length(proportions)
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid prior")
    if(length(prior) != ng) stop("prior is of incorrect length")
    names(prior) <- names(counts) <- lev1
    method <- match.arg(method)
    if(CV && !(method == "moment" || method == "mle"))
        stop(paste("Cannot use leave-one-out CV with method", method))
    group.means <- tapply(x, list(rep(g, p), col(x)), mean)
    f1 <- sqrt(diag(var(x - group.means[g,  ])))
    if(any(f1 < tol))
        stop(paste("variable(s)",
                   paste(format((1:p)[f1 < tol]), collapse = " "),
                   "appear to be constant within groups"))
    # scale columns to unit variance before checking for collinearity
    scaling <- diag(1/f1,,p)
    if(method == "mve") {
        # adjust to "unbiased" scaling of covariance matrix
        cov <- n/(n-ng) * cov.mve((x - group.means[g,  ]) %*% scaling, , FALSE)$cov
        sX <- svd(cov, nu = 0)
        rank <- sum(sX$d > tol^2)
        if(rank < p) warning("variables are collinear")
        scaling <- scaling %*% sX$v[, 1:rank] %*%
            diag(sqrt(1/sX$d[1:rank]),,rank)
    } else if(method == "t") {
        if(nu <= 2) stop("nu must exceed 2")
        w <- rep(1, n)
        repeat {
            w0 <- w
            X <- x - group.means[g, ]
            sX <- svd(sqrt((1 + p/nu)*w/n) * X, nu=0)
            X <- X %*% sX$v %*% diag(1/sX$d,, p)
            w <- 1/(1 + drop(X^2 %*% rep(1, p))/nu)
            print(summary(w))
            group.means <- tapply(w*x, list(rep(g, p), col(x)), sum)/
                rep(tapply(w, g, sum), p)
            if(all(abs(w - w0) < 1e-2)) break
        }
        X <-  sqrt(nu/(nu-2)*(1 + p/nu)/n * w) * (x - group.means[g,  ]) %*% scaling
        X.s <- svd(X, nu = 0)
        rank <- sum(X.s$d > tol)
        if(rank < p) warning("variables are collinear")
        scaling <- scaling %*% X.s$v[, 1:rank] %*% diag(1/X.s$d[1:rank],,rank)
    } else {
        if(method == "moment") fac <- 1/(n-ng) else fac <- 1/n
        X <- sqrt(fac) * (x - group.means[g,  ]) %*% scaling
        X.s <- svd(X, nu = 0)
        rank <- sum(X.s$d > tol)
        if(rank < p) warning("variables are collinear")
        scaling <- scaling %*% X.s$v[, 1:rank] %*% diag(1/X.s$d[1:rank],,rank)
    }
    # now have variables scaled so that W is the identity
    if(CV) {
        x <- x %*% scaling
        dm <- group.means %*% scaling
        K <- if(method == "moment") ng else 0
        dist <- matrix(0, n, ng)
        for(i in 1:ng) {
            dev <- x - matrix(dm[i,  ], n, p, byrow = TRUE)
            dist[, i] <- apply(dev^2, 1, sum)
        }
        ind <- cbind(1:n, g)
        nc <- counts[g]
        cc <- nc/((nc-1)*(n-K))
        dist2 <- dist
        for(i in 1:ng) {
            dev <- x - matrix(dm[i,  ], n, p, byrow = TRUE)
            dev2 <- x - dm[g, ]
            tmp <- apply(dev*dev2, 1, sum)
            dist[, i] <- (n-1-K)/(n-K) * (dist2[, i] +  cc*tmp^2/(1 - cc*dist2[ind]))
        }
        dist[ind] <- dist2[ind] * (n-1-K)/(n-K) * (nc/(nc-1))^2 /
            (1 - cc*dist2[ind])
        dist <- 0.5 * dist - matrix(log(prior), n, ng, byrow=TRUE)
        dist <- exp(-(dist - min(dist, na.rm=TRUE)))
        cl <- factor(apply(dist, 1, which.is.max), levels=seq(along=lev),
                     labels=lev)
        #  convert to posterior probabilities
        posterior <- dist/drop(dist %*% rep(1, length(prior)))
        dimnames(posterior) <- list(rownames(x), lev1)
        return(list(class = cl, posterior = posterior))
    }
    xbar <- apply(prior %*% group.means, 2, sum)
    if(method == "mle") fac <-  1/ng else fac <- 1/(ng - 1)
    X <- sqrt((n * prior)*fac) * scale(group.means, xbar, FALSE) %*% scaling
    X.s <- svd(X, nu = 0)
    rank <- sum(X.s$d > tol * X.s$d[1])
    scaling <- scaling %*% X.s$v[, 1:rank]
    if(is.null(dimnames(x)))
        dimnames(scaling) <- list(NULL, paste("LD", 1:rank, sep = ""))
    else {
        dimnames(scaling) <- list(colnames(x), paste("LD", 1:rank, sep = ""))
        dimnames(group.means)[[2]] <- colnames(x)
    }
    dimnames(group.means)[[1]] <- names(prior)
    structure(list(prior = prior, counts = counts, means = group.means,
                   scaling = scaling, lev = lev, svd = X.s$d[1:rank],
                   N = n, call = match.call()),
              class = "lda")
}

predict.lda <- function(object, newdata, prior = object$prior, dimen,
			method = c("plug-in", "predictive", "debiased"), ...)
{
    if(!inherits(object, "lda")) stop("object not of class lda")
    which.is.max <- function(x) {
        d <- (1:length(x))[x == max(x)]
        if(length(d) > 1) d <- sample(d, 1)
        d
    }
    model.frame.lda <- function(formula)
    {
        oc <- formula$call
        oc$method <- "model.frame"
        oc[[1]] <- as.name("lm")
        eval(oc, sys.frame(sys.parent()))
    }

    if(!is.null(Terms <- object$terms)) {#
        # formula fit
        if(missing(newdata)) newdata <- model.frame.lda(object)
        else newdata <- model.frame(as.formula(delete.response(Terms)),
                                    newdata, na.action=function(x) x)
        x <- model.matrix(delete.response(Terms), newdata)
        xint <- match("(Intercept)", colnames(x), nomatch=0)
        if(xint > 0) x <- x[, -xint, drop=FALSE]
    } else {                            #
        # matrix or data-frame fit
        if(missing(newdata)) {
            if(!is.null(sub <- object$call$subset))
                newdata <- eval(parse(text=paste(deparse(object$call$x),"[",
                                      deparse(sub),",]")), sys.frame(sys.parent()))
            else newdata <- eval(object$call$x, sys.frame(sys.parent()))
            if(!is.null(nas <- object$call$na.action))
                newdata <- eval(call(nas, newdata))
        }
        if(is.null(dim(newdata)))
            dim(newdata) <- c(1, length(newdata))# a row vector
        x <- as.matrix(newdata)		# to cope with dataframes
    }

    if(ncol(x) != ncol(object$means)) stop("wrong number of variables")
    if(length(colnames(x)) > 0 &&
       any(colnames(x) != dimnames(object$means)[[2]]))
        warning("Variable names in newdata do not match those in object")
    ng <- length(prior)
    #   remove overall means to keep distances small
    means <- apply(object$means, 2, mean)
    scaling <- object$scaling
    x <- scale(x, means, FALSE) %*% scaling
    dm <- scale(object$means, means, FALSE) %*% scaling
    method <- match.arg(method)
    if(missing(dimen)) dimen <- length(object$svd)
    else dimen <- min(dimen, length(object$svd))
    N <- object$N
    if(method == "plug-in") {
        dm <- dm[, 1:dimen, drop=FALSE]
        dist <- matrix(0.5 * apply(dm^2, 1, sum) - log(prior), nrow(x),
                       length(prior), byrow = TRUE) - x[, 1:dimen, drop=FALSE] %*% t(dm)
        dist <- exp( -(dist - min(dist, na.rm=TRUE)))
    } else if (method == "debiased") {
        dm <- dm[, 1:dimen, drop=FALSE]
        dist <- matrix(0.5 * apply(dm^2, 1, sum), nrow(x), ng, byrow = TRUE) -
            x[, 1:dimen, drop=FALSE] %*% t(dm)
        dist <- (N - ng - dimen - 1)/(N - ng) * dist -
            matrix(log(prior) - dimen/object$counts , nrow(x), ng, byrow=TRUE)
        dist <- exp( -(dist - min(dist, na.rm=TRUE)))
    } else {                            # predictive
        dist <- matrix(0, nrow = nrow(x), ncol = ng)
        p <- ncol(object$means)
        # adjust to ML estimates of covariances
        X <- x * sqrt(N/(N-ng))
        for(i in 1:ng) {
            nk <- object$counts[i]
            dev <- scale(X, dm[i, ], FALSE)
            dev <- 1 + apply(dev^2, 1, sum) * nk/(N*(nk+1))
            dist[, i] <- prior[i] * (nk/(nk+1))^(p/2) * dev^(-(N - ng + 1)/2)
        }
    }
    cl <- factor(apply(dist, 1, which.is.max), levels=seq(along=object$lev),
                 labels=object$lev)
    posterior <- dist / drop(dist %*% rep(1, ng))
    dimnames(posterior) <- list(rownames(x), names(prior))
    list(class = cl, posterior = posterior, x = x[, 1:dimen, drop=FALSE])
}

print.lda <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }
    cat("\nPrior probabilities of groups:\n")
    print(x$prior, ...)
    cat("\nGroup means:\n")
    print(x$means, ...)
    cat("\nCoefficients of linear discriminants:\n")
    print(x$scaling, ...)
    svd <- x$svd
    names(svd) <- dimnames(x$scaling)[[2]]
    if(length(svd) > 1) {
        cat("\nProportion of trace:\n")
        print(round(svd^2/sum(svd^2), 4), ...)
    }
    invisible(x)
}

plot.lda <- function(x, panel=panel.lda, ..., cex=0.7,
                     dimen, abbrev=FALSE, xlab="LD1", ylab="LD2")
{
    panel.lda <- function(x, y, ...) {
        text(x, y, as.character(g.lda), cex=tcex)
    }
    model.frame.lda <- function(formula)
    {
        oc <- formula$call
        oc$method <- "model.frame"
        oc[[1]] <- as.name("lm")
        eval(oc, sys.frame(sys.parent()))
    }
    if(!is.null(Terms <- x$terms)) {    #
        # formula fit
        data <- model.frame.lda(x)
        X <- model.matrix(delete.response(Terms), data)
        g <- model.extract(data, "response")
        xint <- match("(Intercept)", colnames(X), nomatch=0)
        if(xint > 0) X <- X[, -xint, drop=FALSE]
    } else {                            #
        # matrix or data-frame fit
        xname <- x$call$x
        gname <- x$call[[3]]
        if(!is.null(sub <- x$call$subset)) {
            X <- eval(parse(text=paste(deparse(xname),"[", deparse(sub),",]")),
                      sys.frame(sys.parent()))
            g <- eval(parse(text=paste(deparse(gname),"[", deparse(sub),"]")),
                      sys.frame(sys.parent()))
        } else {
            X <- eval(xname, sys.frame(sys.parent()))
            g <- eval(gname, sys.frame(sys.parent()))
        }
        if(!is.null(nas <- x$call$na.action)) {
            df <- data.frame(g = g, X = X)
            df <- eval(call(nas, df))
            g <- df$g
            X <- df$X
        }
    }
    if(abbrev) levels(g) <- abbreviate(levels(g), abbrev)
    assign("g.lda", g)
    assign("tcex", cex)
    means <- apply(x$means, 2, mean)
    X <- scale(X, means, FALSE) %*% x$scaling
    if(!missing(dimen) && dimen < ncol(X)) X <- X[, 1:dimen, drop=FALSE]
    if(ncol(X) > 2) {
        pairs(X, panel=panel, ...)
    } else if(ncol(X) == 2)  {
        eqscplot(X[, 1:2], xlab=xlab, ylab=ylab, type="n", ...)
        panel(X[, 1], X[, 2], ...)
    } else ldahist(X[,1], g, xlab=xlab, ...)
    invisible(NULL)
}

ldahist <-
function(data, g, nbins = 25, h, x0 = -h/1000, breaks,
	 xlim = range(breaks), ymax = 0, width,
         type = c("histogram", "density", "both"), sep = (type != "density"),
         col = 5,
	 xlab = deparse(substitute(data)), bty = "n", ...)
{
    eval(xlab)
    type <- match.arg(type)
    data <- data[!is.na(data)]
    g <- g[!is.na(data)]
    counts <- table(g)
    groups <- names(counts)[counts > 0]
    if(missing(breaks)) {
        if(missing(h)) h <- diff(pretty(data, nbins))[1]
        first <- floor((min(data) - x0)/h)
        last <- ceiling((max(data) - x0)/h)
        breaks <- x0 + h * c(first:last)
    }
    if(type=="histogram" || type=="both") {
        if(any(diff(breaks) <= 0)) stop("breaks must be strictly increasing")
        if(min(data) < min(breaks) || max(data) > max(breaks))
            stop("breaks do not cover the data")
        est <- vector("list", length(groups))
        for (grp in groups){
            bin <- cut(data[g==grp], breaks, include.lowest = TRUE)
            est1 <- tabulate(bin, length(levels(bin)))
            est1 <- est1/(diff(breaks) * length(data[g==grp]))
            ymax <- max(ymax, est1)
            est[[grp]] <- est1
        }
    }
    if(type=="density" || type == "both"){
        xd <- vector("list", length(groups))
        for (grp in groups){
            if(missing(width)) width <- width.SJ(data[g==grp])
            xd1 <- density(data[g==grp], n=200, width=width,
                           from=xlim[1], to=xlim[2])
            ymax <- max(ymax, xd1$y)
            xd[[grp]] <- xd1
        }
    }
    if(!sep) plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = "",
                  bty = bty)
    else {
        oldpar <- par(mfrow=c(length(groups), 1))
        on.exit(par(oldpar))
    }
    for (grp in groups) {
        if(sep) plot(xlim, c(0, ymax), type = "n",
                     xlab = paste("group", grp), ylab = "", bty = bty)
        if(type=="histogram" || type=="both") {
            n <- length(breaks)
            rect(breaks[-n], 0, breaks[-1], est[[grp]], col = col, ...)
        }
        if(type=="density" || type == "both") lines(xd[[grp]])
    }
    invisible()
}

pairs.lda <- function(x, labels = colnames(x), panel=panel.lda,
                      dimen, abbrev=FALSE, ..., cex=0.7,
                      type=c("std", "trellis"))
{
    panel.lda <- function(x,y, ...) {
        text(x, y, as.character(g.lda), cex=tcex, ...)
    }
    model.frame.lda <- function(formula)
    {
        oc <- formula$call
        oc$method <- "model.frame"
        oc[[1]] <- as.name("lm")
        eval(oc, sys.frame(sys.parent()))
    }
    type <- match.arg(type)
    if(!is.null(Terms <- x$terms)) {    #
        # formula fit
        data <- model.frame.lda(x)
        X <- model.matrix(delete.response(Terms), data)
        g <- model.extract(data, "response")
        xint <- match("(Intercept)", colnames(X), nomatch=0)
        if(xint > 0) X <- X[, -xint, drop=FALSE]
    } else {                            #
        # matrix or data-frame fit
        xname <- x$call$x
        gname <- x$call[[3]]
        if(!is.null(sub <- x$call$subset)) {
            X <- eval(parse(text=paste(deparse(xname),"[", deparse(sub),",]")),
                      sys.frame(sys.parent()))
            g <- eval(parse(text=paste(deparse(gname),"[", deparse(sub),"]")),
                      sys.frame(sys.parent()))
        } else {
            X <- eval(xname, sys.frame(sys.parent()))
            g <- eval(gname, sys.frame(sys.parent()))
        }
        if(!is.null(nas <- x$call$na.action)) {
            df <- data.frame(g = g, X = X)
            df <- eval(call(nas, df))
            g <- df$g
            X <- df$X
        }
    }
    g <- as.factor(g)
    if(abbrev) levels(g) <- abbreviate(levels(g), abbrev)
    assign("g.lda", g)
    assign("tcex", cex)
    means <- apply(x$means, 2, mean)
    X <- scale(X, means, FALSE) %*% x$scaling
    if(!missing(dimen) && dimen < ncol(X)) X <- X[, 1:dimen]
    if(type=="std") pairs.default(X, panel=panel, ...)
    else {
        print(splom(~X, groups = g, panel=panel.superpose,
                    key = list(
                    text=list(levels(g)),
                    points = Rows(trellis.par.get("superpose.symbol"),
                    seq(along=levels(g))),
                    columns = min(5, length(levels(g)))
                    )
                    ))
    }
    invisible(NULL)
}
