# file MASS/dose.p.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
dose.p <- function(obj, cf = 1:2, p = 0.5) {
  eta <- family(obj)$linkfun(p)
  b <- coef(obj)[cf]
  x.p <- (eta - b[1])/b[2]
  names(x.p) <- paste("p = ", format(p), ":", sep = "")
  pd <-  -cbind(1, x.p)/b[2]
  SE <- sqrt(((pd %*% vcov(obj)[cf, cf]) * pd) %*% c(1, 1))
  res <- structure(x.p, SE = SE, p = p)
  class(res) <- "glm.dose"
  res
}

print.glm.dose <- function(x, ...)
{
  M <- cbind(x, attr(x, "SE"))
  dimnames(M) <- list(names(x), c("Dose", "SE"))
  x <- M
  NextMethod("print")
}
