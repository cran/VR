# File MASS/profiles.q copyright (C) 1996 D. M. Bates and W. N. Venables.
#
# port to R by B. D. Ripley copyright (C) 1998
#
#corrections copyright (C) 2000, 3 B. D. Ripley
profile.glm <- function(fitted, which = 1:p, alpha = 0.01,
			maxsteps = 10, del = zmax/5, trace = FALSE, ...)
{
    Pnames <- names(B0 <- coefficients(fitted))
    pv0 <- t(as.matrix(B0))
    p <- length(Pnames)
    if(is.character(which)) which <- match(which, Pnames)
    summ <- summary(fitted)
    std.err <- summ$coefficients[, "Std. Error"]
    mf <- update(fitted, method = "model.frame")
    n <- length(Y <- model.response(mf))
    O <- model.offset(mf)
    if(!length(O)) O <- rep(0, n)
    W <- model.weights(mf)
    if(length(W) == 0) W <- rep(1, n)
    OriginalDeviance <- deviance(fitted)
    DispersionParameter <- summ$dispersion
    X <- model.matrix(fitted)
    fam <- family(fitted)
    switch(fam$family,
           binomial = {
               if(!is.null(dim(Y))) {
                   n <- n/2
                   O <- O[1:n]
                   Y <- Y[, 1]/(W <- drop(Y %*% c(1, 1)))
               }
               zmax <- sqrt(qchisq(1 - alpha/2, p))
               profName <- "z"
           },
           poisson = ,
           "Negative Binomial" = {
               zmax <- sqrt(qchisq(1 - alpha/2, p))
               profName <- "z"
           }
           ,
           gaussian = ,
           quasi = ,
           "inverse.gaussian" = ,
           quasibinomial = ,
           quasipoisson = ,
       {
	   zmax <- sqrt(p * qf(1 - alpha/2, p, n - p))
	   profName <- "tau"
       }
           )
    prof <- vector("list", length=length(which))
    names(prof) <- Pnames[which]
    for(i in which) {
        zi <- 0
        pvi <- pv0
        Xi <- X[,  - i, drop = FALSE]
        pi <- Pnames[i]
        for(sgn in c(-1, 1)) {
            if(trace) cat("\nParameter:", pi, c("down", "up")[(sgn + 1)/2 + 1], "\n")
            step <- 0
            z <- 0
            ## LP is the linear predictor including offset.
            LP <- X %*% fitted$coef + O
            while((step <- step + 1) < maxsteps && abs(z) < zmax) {
                bi <- B0[i] + sgn * step * del * std.err[i]
                o <- O + X[, i] * bi
                ## call to glm.fit.null not needed from 1.4.1 on
                fm <- glm.fit(x = Xi, y = Y, weights = W, etastart = LP,
                              offset = o, family = fam,
                              control = fitted$control)
                LP <- Xi %*% fm$coef + o
                ri <- pv0
                ri[, names(coef(fm))] <- coef(fm)
                ri[, pi] <- bi
                pvi <- rbind(pvi, ri)
                zz <- (fm$deviance - OriginalDeviance)/DispersionParameter
                if(zz > - 1e-3) zz <- max(zz, 0)
                else stop("profiling has found a better solution, so original fit had not converged")
                z <- sgn * sqrt(zz)
                zi <- c(zi, z)
            }
        }
        si <- order(zi)
        prof[[pi]] <- structure(data.frame(zi[si]), names = profName)
        prof[[pi]]$par.vals <- pvi[si, ]
    }
    val <- structure(prof, original.fit = fitted, summary = summ)
    class(val) <- c("profile.glm", "profile")
    val
}

plot.profile <-
  ## R version: non-Trellis-based replacement for plot.profile
  function(x, nseg, ...)
{
    nulls <- sapply(x, is.null)
    if (all(nulls)) return(NULL)
    x <- x[!nulls]
    nm <- names(x)
    nr <- ceiling(sqrt(length(nm)))
    oldpar <- par(mfrow = c(nr, nr))
    on.exit(par(oldpar))
    for(nm in names(x)) {
        tau <- x[[nm]][[1]]
        parval <- x[[nm]][[2]][, nm]
        plot(parval, tau, xlab = nm, ylab = "tau", type="n")
        ## allow for profiling failures
        if(sum(tau == 0) == 1) points(parval[tau == 0], 0, pch = 3)
        splineVals <- spline(parval, tau)
        lines(splineVals$x, splineVals$y)
    }
}

pairs.profile <-
  ## Another plot method for profile objects showing pairwise traces.
  ## Recommended only for diagnostic purposes.
function(x, colours = 2:3, ...)
{
    parvals <- lapply(x, "[[", "par.vals")
    rng <- apply(do.call("rbind", parvals), 2, range, na.rm = TRUE)
    Pnames <- colnames(rng)
    npar <- length(Pnames)
    coefs <- coef(attr(x, "original.fit"))
    form <- paste(as.character(formula(attr(x, "original.fit")))[c(2, 1, 3)],
                  collapse = "")
    oldpar <- par(mar = c(0, 0, 0, 0), mfrow = c(1, 1),
                  oma = c(3, 3, 6, 3), las = 1)
    on.exit(par(oldpar))
    ##
    ## The following dodge ensures that the plot region is square
    ##
    fin <- par("fin")
    dif <- (fin[2] - fin[1])/2
    if(dif > 0) adj <- c(dif, 0, dif, 0)
    else adj <- c(0,  - dif, 0,  - dif)
    par(omi = par("omi") + adj)
    ##
    ##
    cex <- 1 + 1/npar
    frame()
    mtext(form, side = 3, line = 3, cex = 1.5, outer = TRUE)
    del <- 1/npar
    for(i in 1:npar) {
        ci <- npar - i
        pi <- Pnames[i]
        for(j in 1:npar) {
            pj <- Pnames[j]
            par(fig = del * c(j - 1, j, ci, ci + 1))
            if(i == j) {
                par(new=TRUE)
                plot(rng[, pj], rng[, pi], axes = FALSE,
                     xlab = "", ylab = "", type = "n")
                op <- par(usr = c(-1, 1, -1, 1))
                text(0, 0, pi, cex = cex, adj = 0.5)
                par(op)
            } else {
                col <- colours
                if(i < j) col <- col[2:1]
                if(!is.null(parvals[[pj]])) {
                    par(new=TRUE)
                    plot(spline(x <- parvals[[pj]][, pj],
                                y <- parvals[[pj]][, pi]),
                         type = "l", xlim = rng[, pj],
                         ylim = rng[, pi], axes = FALSE,
                         xlab = "", ylab = "", col = col[2])
                    pu <- par("usr")
                    smidge <- 2/100 * (pu[4] - pu[3])
                    segments(x, pmax(pu[3], y - smidge), x,
                             pmin(pu[4], y + smidge))
                } else
                plot(rng[, pj], rng[, pi], axes = FALSE,
                     xlab = "", ylab = "", type = "n")
                if(!is.null(parvals[[pi]])) {
                    lines(x <- parvals[[pi]][, pj], y <- parvals[[pi]][, pi],
                          type = "l", col = col[1])
                    pu <- par("usr")
                    smidge <- 2/100 * (pu[2] - pu[1])
                    segments(pmax(pu[1], x - smidge), y, pmin(pu[2], x + smidge), y)
                }
                points(coefs[pj], coefs[pi], pch = 3, cex = 3)
            }
            if(i == npar) axis(1)
            if(j == 1) axis(2)
            if(i == 1) axis(3)
            if(j == npar) axis(4)
        }
    }
    par(fig = c(0, 1, 0, 1))
    invisible(x)
}
