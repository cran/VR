# file MASS/R/vcov.R
# copyright (C) 1994-2002 W. N. Venables and B. D. Ripley
#
if(R.version$major == 1 && R.version$minor < 6) {
    ## included in R 1.6.0 and later.
vcov <- function(object, ...) UseMethod("vcov")

vcov.nls <- function(object, ...)
{
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}

vcov.glm <- function(object, ...)
{
    so <- summary(object, corr=FALSE, ...)
    so$dispersion * so$cov.unscaled
}

vcov.lm <- function(object, ...)
{
    so <- summary(object, corr=FALSE)
    so$sigma^2 * so$cov.unscaled
}

vcov.coxph <- vcov.survreg <- function (object,...) object$var

vcov.gls <- function (object, ...) object$varBeta

vcov.lme <- function (object, ...) object$varFix
}
