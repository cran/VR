# file MASS/eqscplot.q
# copyright (C) 1994-2005 W. N. Venables and B. D. Ripley
#
eqscplot <- function(x, y, ratio = 1, tol = 0.04, uin, ...)
{
  dots <- list(...); nmdots <- names(dots)
  Call <- match.call()
  Call$ratio <- Call$tol <- Call$uin <- NULL
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
  Call$x <- x; Call$y <- y
  Call$xlab <- if("xlab" %in% nmdots) dots$xlab else xlab0
  Call$ylab <- if("ylab" %in% nmdots) dots$ylab else ylab0
  xlim <- if("xlim" %in% nmdots) dots$xlim else range(x[is.finite(x)])
  ylim <- if("ylim" %in% nmdots) dots$ylim else range(y[is.finite(y)])
  midx <- 0.5 * (xlim[2] + xlim[1])
  xlim <- midx + (1 + tol) * 0.5 * c(-1, 1) * (xlim[2] - xlim[1])
  midy <- 0.5 * (ylim[2] + ylim[1])
  ylim <- midy + (1 + tol) * 0.5 * c(-1, 1) * (ylim[2] - ylim[1])
  oldpin <- par("pin")
  xuin <- oxuin <- oldpin[1]/abs(diff(xlim))
  yuin <- oyuin <- oldpin[2]/abs(diff(ylim))
  if(missing(uin)) {
    if(yuin > xuin*ratio) yuin <- xuin*ratio
    else xuin <- yuin/ratio
  } else {
    if(length(uin) == 1) uin <- uin * c(1, ratio)
    if(any(c(xuin, yuin) < uin)) stop("'uin' is too large to fit plot in")
    xuin <- uin[1]; yuin <- uin[2]
  }
  xlim <- midx + oxuin/xuin * c(-1, 1) * diff(xlim) * 0.5
  ylim <- midy + oyuin/yuin * c(-1, 1) * diff(ylim) * 0.5
  Call$xlim <- xlim
  Call$ylim <- ylim
  Call$xaxs <- Call$yaxs <- "i"
  Call[[1]] <- as.name("plot")
  #plot(x, y, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
  #     xlab = xlab, ylab = ylab, ...)
  eval.parent(Call)
}
