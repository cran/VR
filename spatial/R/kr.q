# file spatial/kr.q copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
.spfmat <- function(x, y, np)
{
    n <- length(x)
    npar <- ((np + 1) * (np + 2))/2
    .C("VR_fmat",
       f = double(n * npar),
       as.double(x),
       as.double(y),
       as.integer(n),
       as.integer(np))$f
}

surf.ls <- function(np, x, y, z)
{
    if(!is.loaded(symbol.C("VR_frset")))
        stop("Compiled code has not been dynamically loaded")
    if (np > 6) stop("np exceeds 6")
    if(is.data.frame(x)) {
        if(any(is.na(match(c("x", "y", "z"), names(x)))))
            stop("`x' does not have columns `x', `y' and `z'")
        if(missing(y)) y <- x$y
        if(missing(z)) z <- x$z
        x <- x$x
    }
    rx <- range(x)
    ry <- range(y)
    .C("VR_frset",
       as.double(rx[1]),
       as.double(rx[2]),
       as.double(ry[1]),
       as.double(ry[2])
       )
    n <- length(x)
    npar <- ((np + 1) * (np + 2))/2
    f <- .spfmat(x, y, np)
    Z <- .C("VR_ls",
            as.double(x),
            as.double(y),
            as.double(z),
            as.integer(n),
            as.integer(np),
            as.integer(npar),
            f = f,
            r = double((npar * (npar + 1))/2),
            beta = double(npar),
            wz = double(n),
            ifail = as.integer(0))
    res <- list(x=x, y=y, z=z, np=np, f=f, r=Z$r, beta=Z$beta,
                wz=Z$wz, rx=rx, ry=ry, call=match.call())
    class(res) <- "trls"
    res
}


surf.gls <- function(np, covmod, x, y, z, nx=1000, ...)
{
    if (np > 6) stop("np exceeds 6")
    if(!is.loaded(symbol.C("VR_frset")))
        stop("Compiled code has not been dynamically loaded")
    if(is.data.frame(x)) {
        if(any(is.na(match(c("x", "y", "z"), names(x)))))
            stop("`x' does not have columns `x', `y' and `z'")
        if(missing(y)) y <- x$y
        if(missing(z)) z <- x$z
        x <- x$x
    }
    rx <- range(x)
    ry <- range(y)
    .C("VR_frset", as.double(rx[1]), as.double(rx[2]),
       as.double(ry[1]), as.double(ry[2]))
    covmod <- covmod
    arguments <- list(...)
    if (length(arguments)) {
        onames <- names(formals(covmod))
        pm <- pmatch(names(arguments), onames, nomatch = 0)
        if (any(pm == 0)) warning(paste("some of ... do not match"))
        names(arguments[pm > 0]) <- onames[pm]
        oargs <- formals(covmod)
        oargs[pm] <- arguments[pm > 0]
        formals(covmod) <- oargs
    }
    mm <- 1.5*sqrt((rx[2]-rx[1])^2 + (ry[2]-ry[1])^2)
    alph <- c(mm/nx, covmod(seq(0, mm, mm/nx)))
    .C("VR_alset", as.double(alph), as.integer(length(alph)))
    n <- length(x)
    npar <- ((np + 1) * (np + 2))/2
    f <- .spfmat(x, y, np)
    Z <- .C("VR_gls",
            as.double(x),
            as.double(y),
            as.double(z),
            as.integer(n),
            as.integer(np),
            as.integer(npar),
            as.double(f),
            l = double((n * (n + 1))/2),
            r = double((npar * (npar + 1))/2),
            beta = double(npar),
            wz = double(n),
            yy = double(n),
            W = double(n),
            ifail = as.integer(0),
            l1f = double(n * npar)
            )
    if(Z$ifail > 0) stop("Rank failure in Choleski decomposition")
    if(nx > 1000) alph <- alph[1]
    res <- list(x=x, y=y, z=z, np=np, f=f, alph=alph, l=Z$l, r=Z$r,
                beta=Z$beta, wz=Z$wz, yy=Z$yy, W=Z$W, l1f=Z$l1f, rx=rx, ry=ry,
                covmod=covmod, call=match.call())
    class(res) <- c("trgls", "trls")
    res
}

.trval <- function(obj, x, y)
{
    n <- length(x)
    .C("VR_valn",
       z = double(n),
       as.double(x),
       as.double(y),
       as.integer(n),
       as.double(obj$beta),
       as.integer(obj$np))$z
}

predict.trls <- function (object, x, y, ...)
{
    if (!inherits(object, "trls"))
        stop("object not a fitted trend surface")
    n <- length(x)
    if (length(y) != n) stop("x and y differ in length")
    .C("VR_frset",
       as.double(object$rx[1]),
       as.double(object$rx[2]),
       as.double(object$ry[1]),
       as.double(object$ry[2])
       )
    invisible(.trval(object, x, y))
}

trmat <- function (obj, xl, xu, yl, yu, n)
{
    if (!inherits(obj, "trls"))
        stop("object not a fitted trend surface")
    dx <- (xu - xl)/n
    dy <- (yu - yl)/n
    x <- seq(xl, xu, dx)
    y <- seq(yl, yu, dy)
    z <- matrix(nrow = length(x), ncol = length(y))
    for (i in seq(along = y))
        z[, i] <- predict(obj, x, rep(y[i], length(x)))
    invisible(list(x = x, y = y, z = z))
}

if(0){
trmat <- function(obj, xl, xu, yl, yu, n)
{
    if(!inherits(obj, "trls")) stop("object not a fitted trend surface")
    .C("VR_frset",
       as.double(obj$rx[1]),
       as.double(obj$rx[2]),
       as.double(obj$ry[1]),
       as.double(obj$ry[2])
       )
    dx <- (xu - xl)/n
    dy <- (yu - yl)/n
    x <- seq(xl, xu, dx)
    y <- seq(yl, yu, dy)
    z <- matrix(nrow = length(x), ncol = length(y))
    for(i in seq(along = y))
        z[, i] <- .trval(obj, x, rep(y[i], length(x)))
    invisible(list(x = x, y = y, z = z))
}
}

prmat <- function(obj, xl, xu, yl, yu, n)
{
    predval <- function(obj, xp, yp)
    {
        npt <- length(xp)
        .C("VR_krpred",
           z = double(npt),
           as.double(xp),
           as.double(yp),
           as.double(obj$x),
           as.double(obj$y),
           as.integer(npt),
           as.integer(length(obj$x)),
           as.double(obj$yy)
           )$z
    }

    if(!inherits(obj, "trgls")) stop("object not from kriging")
    if(n > 999) stop("n is too large")
    .C("VR_frset",
       as.double(obj$rx[1]),
       as.double(obj$rx[2]),
       as.double(obj$ry[1]),
       as.double(obj$ry[2])
       )
    alph <- obj$alph
    if(length(alph) <= 1) {
        mm <- 1.5*sqrt((obj$rx[2]-obj$rx[1])^2 + (obj$ry[2]-obj$ry[1])^2)
        alph <- c(alph[1], obj$covmod(seq(0, mm, alph[1])))
    }
    .C("VR_alset",
       as.double(alph),
       as.integer(length(alph))
       )
    dx <- (xu - xl)/n
    dy <- (yu - yl)/n
    xs <- seq(xl, xu, dx)
    ys <- seq(yl, yu, dy)
    z <- matrix(nrow = length(xs), ncol = length(ys))
    for(i in seq(along = ys))
        z[, i] <- .trval(obj, xs, rep(ys[i], length(xs))) +
            predval(obj, xs, rep(ys[i], length(xs)))
    invisible(list(x = xs, y = ys, z = z))
}


semat <- function(obj, xl, xu, yl, yu, n, se)
{
    seval <- function(obj, xp, yp)
    {
        npt <- length(xp)
        np <- obj$np
        npar <- ((np + 1) * (np + 2))/2
        .C("VR_prvar",
           z = double(npt),
           as.double(xp),
           as.double(yp),
           as.integer(npt),
           as.double(obj$x),
           as.double(obj$y),
           as.double(obj$l),
           as.double(obj$r),
           as.integer(length(obj$x)),
           as.integer(np),
           as.integer(npar),
           as.double(obj$l1f)
           )$z
    }

    if(!inherits(obj, "trgls")) stop("object not from kriging")
    if(n > 999) stop("n is too large")
    .C("VR_frset",
       as.double(obj$rx[1]),
       as.double(obj$rx[2]),
       as.double(obj$ry[1]),
       as.double(obj$ry[2])
       )
    alph <- obj$alph
    if(length(alph) <= 1) {
        mm <- 1.5*sqrt((obj$rx[2]-obj$rx[1])^2 + (obj$ry[2]-obj$ry[1])^2)
        alph <- c(alph[1], obj$covmod(seq(0, mm, alph[1])))
  }
    .C("VR_alset", as.double(alph), as.integer(length(alph)))
    dx <- (xu - xl)/n
    dy <- (yu - yl)/n
    xs <- seq(xl, xu, dx)
    ys <- seq(yl, yu, dy)
    z <- matrix(nrow = length(xs), ncol = length(ys))
    np <- obj$np
    npar <- ((np + 1) * (np + 2))/2
    if(missing(se))
        se <- sqrt(sum(obj$W^2)/(length(obj$x) - npar))
    for(i in seq(along = ys))
        z[, i] <- se * sqrt(seval(obj, xs, rep(ys[i], length(xs))))
    invisible(list(x = xs, y = ys, z = z))
}

correlogram <- function(krig, nint, plotit=TRUE, ...)
{
  if(!is.loaded(symbol.C("VR_correlogram")))
    stop("Compiled code has not been dynamically loaded")
  z <- .C("VR_correlogram",
          xp = double(nint),
          yp = double(nint),
          nint = as.integer(nint),
          as.double(krig$x),
          as.double(krig$y),
          if(krig$np > 0) as.double(krig$wz) else as.double(krig$z),
          as.integer(length(krig$x)),
          cnt = integer(nint)
          )
  xp <- z$xp[1:z$nint]
  yp <- z$yp[1:z$nint]
  z <- list(x = xp, y = yp, cnt = z$cnt[1:z$nint])
  if(plotit)
    if(exists(".Device")) {
      plot(xp, yp, type = "p", ylim = c(-1, 1), ...)
      abline(0, 0)
      invisible(z)
    }
    else {
      warning("Device not active")
      return(z)
    }
  else z
}

variogram <- function(krig, nint, plotit=TRUE, ...)
{
  if(!is.loaded(symbol.C("VR_variogram")))
    stop("Compiled code has not been dynamically loaded")
  z <- .C("VR_variogram",
          xp = double(nint),
          yp = double(nint),
          nint = as.integer(nint),
          as.double(krig$x),
          as.double(krig$y),
          if(krig$np > 0) as.double(krig$wz) else as.double(krig$z),
          as.integer(length(krig$x)),
          cnt = integer(nint)
          )
  xp <- z$xp[1:z$nint]
  yp <- z$yp[1:z$nint]
  if(xp[1] > 0) {xp <- c(0, xp); yp <- c(0, yp)}
  z <- list(x = xp, y = yp, cnt = z$cnt[1:z$nint])
  if(plotit)
    if(exists(".Device")) {
      plot(xp, yp, type = "p", ...)
      invisible(z)
    }
    else {
      warning("Device not active")
      return(z)
    }
  else z
}

expcov <- function(r, d, alpha=0, se=1)
{
    se^2*(alpha*(r < d/10000) + (1-alpha)*exp(-r/d))
}

gaucov <- function(r, d, alpha=0, se=1)
{
    se^2*(alpha*(r < d/10000) + (1-alpha)*exp(-(r/d)^2))
}

sphercov <- function(r, d, alpha=0, se=1, D=2)
{
    r <- r/d
    if(D == 2) {
        t <- 1 - (2/pi)*(r*sqrt(1-r^2) + asin(r))
    } else {
        t <- 1 - 1.5*r + r^3/2
    }
    se^2*(alpha*(r < 1/10000) + (1-alpha)*ifelse(r < 1, t, 0))
}

# Method despatch functions for trend surface trls class objects
#
residuals.trls <- function (object, ...)
{
    if (!inherits(object, "trls"))
            stop("object not a fitted trend surface")
    object$wz
}

fitted.trls <- function (object, ...)
{
    if (!inherits(object, "trls"))
        stop("object not a fitted trend surface")
    object$z - residuals(object)
}

deviance.trls <- function (object, ...)
{
    if (!inherits(object, "trls"))
        stop("object not a fitted trend surface")
    sum(residuals(object)^2)
}

df.residual.trls <- function (object, ...)
{
    if (!inherits(object, "trls"))
        stop("object not a fitted trend surface")
    length(object$z) - length(object$beta)
}

extractAIC.trls <- function (fit, scale, k = 2, ...)
{
    if (!inherits(fit, "trls"))
        stop("object not a fitted trend surface")
    n <- length(fit$z)
    edf <- df.residual.trls(fit)
    RSS <- deviance(fit)
    dev <- n * log(RSS/n)
    c(edf, dev + k * edf)
}

#
# Anova output to match Burrough & McDonnell (1998) Principals
# of Geographical Information Systems (Oxford University Press)
# book tabulation
#
anova.trls <- function (object, ...)
{
    if (length(list(object, ...)) > 1)
        return(anovalist.trls(object, ...))
    if (!inherits(object, "trls"))
        stop("object not a fitted trend surface")
    rss <- deviance(object)
    rdf <- df.residual.trls(object)
    n <- length(object$z)
    edf <- n - rdf - 1
    tss <- var(object$z) * (n - 1)
    ess <- tss - rss
    ems <- ess/edf
    rms <- rss/rdf
    f <- ems/rms
    p <- 1 - pf(f, edf, rdf)
    table <- data.frame(format(c(ess, rss, tss)),
                        format(c(edf, rdf, edf + rdf)),
                        c(format(c(ems, rms)), ""),
                        c(format(f), "", ""),
                        c(format.pval(p), "", ""))
    dimnames(table) <-
        list(c("Regression", "Deviation", "Total"),
             c("Sum Sq", "Df", "Mean Sq", "F value", "Pr(>F)"))
    cat("Analysis of Variance Table\n", "Model: ")
    cat(deparse(object$call), "\n", sep="")
    table
}

anovalist.trls <- function (object, ...)
{
    objs <- list(object, ...)
    nmodels <- length(objs)
    for (i in 1:nmodels) {
        if (!inherits(objs[[i]], "trls"))
            stop("object not a fitted trend surface")
    }
    if (nmodels == 1)
        return(anova.trls(object))
    models <- as.character(lapply(objs, function(x) x$call))
    df.r <- unlist(lapply(objs, df.residual.trls))
    ss.r <- unlist(lapply(objs, deviance.trls))
    df <- c(NA, -diff(df.r))
    ss <- c(NA, -diff(ss.r))
    ms <- ss/df
    f <- p <- rep(NA, nmodels)
    for (i in 2:nmodels) {
        if (df[i] > 0) {
            f[i] <- ms[i]/(ss.r[i]/df.r[i])
            p[i] <- 1 - pf(f[i], df[i], df.r[i])
        } else if (df[i] < 0) {
            f[i] <- ms[i]/(ss.r[i - 1]/df.r[i - 1])
            p[i] <- 1 - pf(f[i], -df[i], df.r[i - 1])
        } else {
            ss[i] <- 0
        }
    }
    table <- data.frame(df.r, ss.r, df, ss, f, p)
    dimnames(table) <-
        list(1:nmodels, c("Res.Df", "Res.Sum Sq",
                          "Df", "Sum Sq", "F value", "Pr(>F)"))
    title <- "Analysis of Variance Table\n"
    topnote <- paste("Model ", format(1:nmodels), ": ", models,
                     sep = "", collapse = "\n")
    sss <- getOption("show.signif.stars")
    if (sss) options(show.signif.stars = FALSE)
    print(structure(table, heading = c(title, topnote),
                    class = c("anova", "data.frame")))
    if (sss) options(show.signif.stars = TRUE)
    invisible(structure(table, heading = c(title, topnote),
                        class = c("anova", "data.frame")))
}

#
# and a basic summary method, avoiding the coefficient values
#
summary.trls <-
    function (object, digits = max(3, getOption("digits") - 3), ...)
{
    if (!inherits(object, "trls"))
        stop("object not a fitted trend surface")
    print(anova.trls(object))
    rdf <- df.residual.trls(object)
    n <- length(object$z)
    edf <- n - rdf - 1
    rss <- deviance(object)
    tss <- var(object$z) * (n - 1)
    ess <- tss - rss
    ems <- ess/edf
    rsquared <- ess/tss
    adj.rsquared <- 1 - (1 - rsquared) * ((n - 1)/rdf)
    cat("Multiple R-Squared:", format(rsquared, digits = digits))
    cat(",\tAdjusted R-squared:", format(adj.rsquared, digits = digits),
        "\n")
    AIC <- extractAIC(object)
    cat("AIC: (df = ", AIC[1], ") ", AIC[2], "\n", sep = "")
    cat("Fitted:\n")
    if (rdf > 5) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- structure(quantile(fitted.trls(object)), names = nam)
        print(rq, digits = digits)
    } else {
        print(fitted(object), digits = digits)
    }
    cat("Residuals:\n")
    if (rdf > 5) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- structure(quantile(residuals(object)), names = nam)
        print(rq, digits = digits)
    } else {
        print(residuals(object), digits = digits)
    }
}
#
# Influence measures, rather fewer than lm.influence
#
trls.influence <- function (object)
{
    if (!inherits(object, "trls"))
        stop("object not a fitted trend surface")
    nr <- length(object$z)
    nc <- length(object$beta)
    X <- matrix(object$f, nrow = nr, ncol = nc)
    hii <- hat(X, FALSE)
    s <- sqrt(deviance(object)/df.residual.trls(object))
    r <- residuals(object)
    stresid <- r/(s * sqrt(1 - hii))
    Di <- ((stresid^2) * hii)/(nc * (1 - hii))
    invisible(list(r = r, hii = hii, stresid = stresid, Di = Di))
}

plot.trls <-
    function (x, border = "red", col = NA, pch = 4, cex = 0.6,
              add = FALSE, div = 8, ...)
{
    if (!inherits(x, "trls"))
        stop("`x' not a fitted trend surface")
    infl <- trls.influence(x)
    dx <- diff(range(x$x))
    dy <- diff(range(x$y))
    dxy <- (dx + dy)/2
    mDi <- max(infl$Di)
    sc <- (mDi * dxy)/div
    if (!add)
        plot(x$x, x$y, type = "n", xlab = "", ylab = "")
    symbols(x$x, x$y, circles=sc * infl$Di, add=TRUE,
            fg=border,
            inches=FALSE)
    points(x$x, x$y, pch = pch)
}

