# file MASS/qda.q
# copyright (C) 1994-2003 W. N. Venables and B. D. Ripley
#
qda <- function(x, ...) UseMethod("qda")

qda.formula <- function(formula, data = NULL, ...,
			subset, na.action = na.fail)
{
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval.parent(m$data)))
        m$data <- as.data.frame(data)
    m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    grouping <- model.extract(m, "response")
    x <- model.matrix(Terms, m)
    xvars <- as.character(attr(Terms, "variables"))[-1]
    if ((yvar <- attr(Terms, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
        xlev <- lapply(m[xvars], levels)
        xlev[!sapply(xlev, is.null)]
    }
    xint <- match("(Intercept)", colnames(x), nomatch=0)
    if(xint > 0) x <- x[, -xint, drop=FALSE]
    res <- qda.default(x, grouping, ...)
    res$terms <- Terms
    res$call <- match.call()
    res$contrasts <- attr(x, "contrasts")
    res$xlevels <- xlev
    attr(res, "na.message") <- attr(m, "na.message")
    if(!is.null(attr(m, "na.action"))) res$na.action <- attr(m, "na.action")
    res
}

qda.data.frame <- function(x, ...)
{
    res <- qda.matrix(structure(data.matrix(x), class="matrix"), ...)
    res$call <- match.call()
    res
}


qda.matrix <- function(x, grouping, ...,
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
    res <- NextMethod("qda")
    res$call <- match.call()
    res
}

qda.default <-
  function(x, grouping, prior = proportions,
           method = c("moment", "mle", "mve", "t"),
           CV = FALSE, nu = 5, ...)
{
    if(is.null(dim(x))) stop("x is not a matrix")
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if(n != length(grouping)) stop("nrow(x) and length(grouping) are different")
    g <- as.factor(grouping)
    lev <- levels(g)
    counts <- as.vector(table(g))
    names(counts) <- lev
    if(any(counts < p+1)) stop("some group is too small for qda")
    proportions <- counts/length(g)
    ng <- length(proportions)
# allow for supplied prior
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid prior")
    if(length(prior) != ng) stop("prior is of incorrect length")
    names(prior) <- lev
# means by group (rows) and variable (columns)
    group.means <- tapply(x, list(rep(g, ncol(x)), col(x)), mean)
    scaling <- array(dim=c(p,p,ng))
    ldet <- numeric(ng)
    method <- match.arg(method)
    if(CV && !(method == "moment" || method == "mle"))
        stop(paste("Cannot use leave-one-out CV with method", method))
    for (i in 1:ng){
        if(method == "mve") {
            cX <- cov.mve(x[unclass(g) == i, ])
            group.means[i,] <- cX$center
            sX <- svd(cX$cov, nu=0)
            scaling[, , i] <- sX$v %*% diag(sqrt(1/sX$d),,p)
            ldet[i] <- sum(log(sX$d))
        } else if(method == "t") {
            if(nu <= 2) stop("nu must exceed 2")
            m <- counts[i]
            X <- x[unclass(g) == i, ]
            w <- rep(1, m)
            repeat {
                w0 <- w
                W <- scale(X, center=group.means[i, ], scale=FALSE)
                sX <- svd(sqrt((1 + p/nu) * w/m) * W, nu=0)
                W <- W %*% sX$v %*% diag(1/sX$d,, p)
                w <- 1/(1 + drop(W^2 %*% rep(1, p))/nu)
                #         print(summary(w))
                group.means[i,] <- colSums(w*X)/sum(w)
                if(all(abs(w - w0) < 1e-2)) break
            }
            qx <- qr(sqrt(w)*scale(X, center=group.means[i, ], scale=FALSE))
            if(qx$rank < p) stop(paste("Rank deficiency in group", lev[i]))
            qx <- qx$qr* sqrt((1 + p/nu)/m)
            scaling[, , i] <- backsolve(qx[1:p,  ], diag(p))
            ldet[i] <- 2*sum(log(abs(diag(qx))))
        } else {
            if(method == "moment") nk <- counts[i] - 1 else nk <- counts[i]
            X <- scale(x[unclass(g) == i, ], center=group.means[i, ], scale=FALSE)/sqrt(nk)
            qx <- qr(X)
            if(qx$rank < p) stop(paste("Rank deficiency in group", lev[i]))
            qx <- qx$qr
            scaling[, , i] <- backsolve(qx[1:p, ], diag(p))
            ldet[i] <- 2*sum(log(abs(diag(qx))))
        }
    }
    if(CV) {
        NG <- if(method == "mle") 0 else 1
        dist <- matrix(0, n, ng)
        Ldet <- matrix(0, n, ng)
        for(i in 1:ng) {
            dev <- ((x - matrix(group.means[i,  ], nrow(x),
                                p, byrow = TRUE)) %*% scaling[,,i])
            dist[, i] <- rowSums(dev^2)
            Ldet[, i] <- ldet[i]
        }
        nc <- counts[g]
        ind <- cbind(1:n, g)
        fac <- 1 - nc/(nc-1)/(nc-NG) * dist[ind]
        fac[] <- pmax(fac, 1e-10)  # possibly degenerate dsn
        Ldet[ind] <- log(fac) + p * log((nc-NG)/(nc-1-NG)) + Ldet[ind]
        dist[ind] <- dist[ind] * (nc^2/(nc-1)^2) * (nc-1-NG)/(nc-NG) / fac
        dist <- 0.5 * dist + 0.5 * Ldet - matrix(log(prior), n, ng, byrow=TRUE)
        dist <- exp(-(dist - min(dist, na.rm=TRUE)))
        posterior <- dist/drop(dist %*% rep(1, length(prior)))
        cl <- factor(max.col(posterior), levels=seq(along=lev), labels=lev)
        dimnames(posterior) <- list(rownames(x), lev)
        return(list(class = cl, posterior = posterior))
    }
    if(is.null(dimnames(x)))
        dimnames(scaling) <- list(NULL, as.character(1:p), lev)
    else {
        dimnames(scaling) <- list(colnames(x), as.character(1:p), lev)
        dimnames(group.means)[[2]] <- colnames(x)
    }
    res <- list(prior = prior, counts = counts, means = group.means,
                scaling = scaling, ldet = ldet, lev = lev, N = n,
                call = match.call())
    class(res) <- "qda"
    res
}

predict.qda <- function(object, newdata, prior = object$prior,
			method = c("plug-in", "predictive", "debiased",
                          "looCV"), ...)
{
    if(!inherits(object, "qda")) stop("object not of class qda")
    method <- match.arg(method)
    if(method == "looCV" && !missing(newdata))
        stop("Cannot have leave-one-out CV with newdata")
    if(is.null(mt <- object$call$method)) mt <- "moment"
    if(method == "looCV" && !(mt == "moment" || mt == "mle"))
        stop(paste("Cannot use leave-one-out CV with method", mt))
    if(!is.null(Terms <- object$terms)) {
    # formula fit
        if(missing(newdata)) newdata <- model.frame(object)
        else {
            newdata <- model.frame(as.formula(delete.response(Terms)),
                                   newdata, na.action=function(x) x,
                                   xlev = object$xlevels)
        }
        x <- model.matrix(delete.response(Terms), newdata,
                          contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(x), nomatch=0)
        if(xint > 0) x <- x[, -xint, drop=FALSE]
        if(method == "looCV") g <- model.extract(newdata, "response")
    } else { #
    # matrix or data-frame fit
        if(missing(newdata)) {
            if(!is.null(sub <- object$call$subset)) {
                newdata <-
                    eval.parent(parse(text=paste(deparse(object$call$x),
                                      "[", deparse(sub),",]")))
                g <- eval.parent(parse(text=paste(deparse(object$call[[3]]),
                                       "[", deparse(sub),"]")))
            } else {
                newdata <- eval.parent(object$call$x)
                g <- eval.parent(object$call[[3]])
            }
            if(!is.null(nas <- object$call$na.action)) {
                df <- data.frame(g = g, X = newdata)
                df <- eval(call(nas, df))
                g <- df$g
                newdata <- df$X
            }
            g <- as.factor(g)
        }
        if(is.null(dim(newdata)))
            dim(newdata) <- c(1, length(newdata))  # a row vector
        x <- as.matrix(newdata)		# to cope with dataframes
    }
    p <- ncol(object$means)
    if(ncol(x) != p) stop("wrong number of variables")
    if(length(colnames(x)) > 0 &&
       any(colnames(x) != dimnames(object$means)[[2]]))
        warning("Variable names in newdata do not match those in object")
    ngroup <- length(object$prior)
    dist <- matrix(0, nrow = nrow(x), ncol = ngroup)
    if(method == "plug-in") {
        for(i in 1:ngroup) {
            dev <- ((x - matrix(object$means[i,  ], nrow(x),
                                ncol(x), byrow = TRUE)) %*% object$scaling[,,i])
            dist[, i] <- 0.5 * rowSums(dev^2) + 0.5 * object$ldet[i] - log(object$prior[i])
        }
#        dist <- exp( -(dist - min(dist, na.rm=T)))
        dist <- exp( -(dist - apply(dist, 1, min, na.rm=TRUE)))
    } else if(method == "looCV") {
        n <- nrow(x)
        NG <- 1
        if(mt == "mle") NG <- 0
        ldet <- matrix(0, n, ngroup)
        for(i in 1:ngroup) {
            dev <- ((x - matrix(object$means[i,  ], nrow(x), p, byrow = TRUE))
                    %*% object$scaling[,,i])
            dist[, i] <- rowSums(dev^2)
            ldet[, i] <- object$ldet[i]
        }
        nc <- object$counts[g]
        ind <- cbind(1:n, g)
        fac <- 1 - nc/(nc-1)/(nc-NG) * dist[ind]
        fac[] <- pmax(fac, 1e-10)  # possibly degenerate dsn
        ldet[ind] <- log(fac) + p * log((nc-NG)/(nc-1-NG)) + ldet[ind]
        dist[ind] <- dist[ind] * (nc^2/(nc-1)^2) * (nc-1-NG)/(nc-NG) / fac
        dist <- 0.5 * dist + 0.5 * ldet -
            matrix(log(object$prior), n, ngroup, byrow=TRUE)
        dist <- exp( -(dist - apply(dist, 1, min, na.rm=TRUE)))
    } else if(method == "debiased") {
        for(i in 1:ngroup) {
            nk <- object$counts[i]
            Bm <- p * log((nk-1)/2) - sum(digamma(0.5 * (nk - 1:ngroup)))
            dev <- ((x - matrix(object$means[i,  ], nrow = nrow(x),
                                ncol = ncol(x), byrow = TRUE)) %*% object$scaling[,,i])
            dist[, i] <- 0.5 * (1 - (p-1)/(nk-1)) * rowSums(dev^2) +
                0.5 * object$ldet[i] - log(object$prior[i]) + 0.5 * Bm - p/(2*nk)
        }
        dist <- exp( -(dist - apply(dist, 1, min, na.rm=TRUE)))
    } else {
        N <- object$N
        for(i in 1:ngroup) {
            nk <- object$counts[i]
            dev <- ((x - matrix(object$means[i,  ], nrow = nrow(x),
                                ncol = ncol(x), byrow = TRUE))
                    %*% object$scaling[,,i])
            dev <- 1 + rowSums(dev^2)/(nk+1)
            dist[, i] <- object$prior[i] * exp(-object$ldet[i]/2) *
                dev^(-nk/2) * (1 + nk)^(-p/2)
        }
    }
    posterior <- dist/drop(dist %*% rep(1, length(object$prior)))
    cl <- factor(max.col(posterior), levels=seq(along=object$lev),
                 labels=object$lev)
    dimnames(posterior) <- list(rownames(x), object$lev)
    list(class = cl, posterior = posterior)
}

print.qda <- function(x, ...)
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
    invisible(x)
}

model.frame.qda <-  model.frame.lda
