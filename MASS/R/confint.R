# file MASS/confint.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
confint <- function(object, parm, level = 0.95, ...) UseMethod("confint")

confint.glm <- function(object, parm, level = 0.95, trace = FALSE, ...)
{
    pnames <- names(coef(object))
    if(missing(parm)) parm <- seq(along=pnames)
    else if(is.character(parm))  parm <- match(parm, pnames, nomatch = 0)
    cat("Waiting for profiling to be done...\n")
    object <- profile(object, which = parm, alpha = (1. - level)/4.,
                      trace = trace)
    confint(object, parm=parm, level=level, trace=trace, ...)
}

confint.profile.glm <-
  function(object, parm = seq(along=pnames), level = 0.95, ...)
{
    of <- attr(object, "original.fit")
    pnames <- names(coef(of))
    if(is.character(parm))  parm <- match(parm, pnames, nomatch = 0)
    a <- (1-level)/2
    a <- c(a, 1-a)
    pct <- paste(round(100*a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2),
                dimnames = list(pnames[parm], pct))
    cutoff <- qnorm(a)
    for(pm in parm) {
        pro <- object[[pm]]
        if(length(pnames) > 1)
            sp <- spline(x = pro[, "par.vals"][, pm], y = pro$z)
        else sp <- spline(x = pro[, "par.vals"], y = pro$z)
        ci[pnames[pm], ] <- approx(sp$y, sp$x, xout = cutoff)$y
    }
    drop(ci)
}

confint.nls <-
  function(object, parm = seq(along=pnames), level = 0.95, ...)
{
  pnames <- names(coef(object))
  if(is.character(parm))  parm <- match(parm, pnames, nomatch = 0)
  cat("Waiting for profiling to be done...\n")
  object <- profile(object, which = parm, alphamax = (1. - level)/4.)
  confint(object, parm=parm, level=level, ...)
}

confint.profile.nls <-
  function(object, parm = seq(along=pnames), level = 0.95, ...)
{
  of <- attr(object, "original.fit")
  pnames <- names(coef(of))
  if(is.character(parm))  parm <- match(parm, pnames, nomatch = 0)
  n <- length(fitted(of)) - length(of$m$getPars())
  a <- (1-level)/2
  a <- c(a, 1-a)
  pct <- paste(round(100*a, 1), "%", sep = "")
  ci <- array(NA, dim = c(length(parm), 2),
              dimnames = list(pnames[parm], pct))
  cutoff <- qt(a, n)
  for(pm in parm) {
    pro <- object[[pm]]
    if(length(pnames) > 1)
        sp <- spline(x = pro[, "par.vals"][, pm], y = pro$tau)
    else sp <- spline(x = pro[, "par.vals"], y = pro$tau)
    ci[pnames[pm], ] <- approx(sp$y, sp$x, xout = cutoff)$y
  }
  drop(ci)
}




