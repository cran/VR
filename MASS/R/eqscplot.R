# file MASS/eqscplot.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
eqscplot <- function(x, y, tol = 0.04, xlim = range(x), ylim = range(y),
		     xlab, ylab,
		     ...)
{
  if(is.matrix(x)) {
    y <- x[, 2]
    x <- x[, 1]
    if(!is.null(dn <- colnames(x))) {
      xlab0 <- dn[1]
      ylab0 <- dn[2]
    } else {
      xlab0 <- ""
      ylab0 <- ""
    }
  } else if(is.list(x)) {
    y <- x$y
    x <- x$x
    xlab0 <- "x"; ylab0 <- "y"
  } else {
    xlab0 <- deparse(substitute(x))
    ylab0 <- deparse(substitute(y))
  }
  if(missing(xlab)) xlab <- xlab0
  if(missing(ylab)) ylab <- ylab0
  oldpin <- par("pin")
  midx <- 0.5 * (xlim[2] + xlim[1])
  xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
  midy <- 0.5 * (ylim[2] + ylim[1])
  ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
  xr <- oldpin[1]/(xlim[2] - xlim[1])
  yr <- oldpin[2]/(ylim[2] - ylim[1])
  if (yr > xr) {
    ylim <- midy + yr * c(-1, 1) * (ylim[2] - ylim[1])/(2*xr)
  } else {
    xlim <- midx + xr * c(-1, 1) * (xlim[2] - xlim[1])/(2*yr)
  }
  plot(x, y, xlim=xlim, ylim=ylim, xaxs="i", yaxs="i",
       xlab=xlab, ylab=ylab, ...)
}
