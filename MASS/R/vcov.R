# file MASS/vcov.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
vcov <- function(object, ...) UseMethod("vcov")

vcov.nls <- function(object, ...)
{
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}

vcov.glm <- function(object, ...)
{
    so <- summary(object, corr=FALSE)
    so$dispersion * so$cov.unscaled
}

vcov.lm <- function(object, ...)
{
    so <- summary(object, corr=FALSE)
    so$sigma^2 * so$cov.unscaled
}

