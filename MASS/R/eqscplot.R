# file MASS/eqscplot.q
# copyright (C) 1994-2001 W. N. Venables and B. D. Ripley
#
eqscplot <- function(x, y, ratio = 1, tol = 0.04, uin,
                     xlim = range(x[is.finite(x)]),
                     ylim = range(y[is.finite(y)]),
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
  midx <- 0.5 * (xlim[2] + xlim[1])
  xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
  midy <- 0.5 * (ylim[2] + ylim[1])
  ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
  oldpin <- par("pin")
  xuin <- oxuin <- oldpin[1]/diff(xlim)
  yuin <- oyuin <- oldpin[2]/diff(ylim)
  if(missing(uin)) {
    if(yuin > xuin*ratio) yuin <- xuin*ratio
    else xuin <- yuin/ratio
  } else {
    if(length(uin) == 1) uin <- uin * c(1, ratio)
    if(any(c(xuin, yuin) < uin)) stop("uin is too large to fit plot in")
    xuin <- uin[1]; yuin <- uin[2]
  }
  xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
  ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5
  plot(x, y, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
       xlab = xlab, ylab = ylab, ...)
}
