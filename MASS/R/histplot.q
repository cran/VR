# file MASS/histplot.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
"histplot" <-
function(formula, data = sys.parent(1),
	 prepanel = prepanel.histplot, panel = panel.histplot,
	 nbins = 5, h , x0 = -h/1000, breaks, prob = TRUE,
	 aspect = "fill",
	 xlab = deparse(do.formula.trellis(formula)$expr[[1]]),
	 ylab = if(prob) "Density" else "Counts",
	 groups = NULL, ..., subset = TRUE)
{
  .NotYetImplemented()
  if(mode(formula) != "call")
    formula <- formula.default(paste("~", deparse(substitute(formula))))
  Z <- do.formula.trellis(formula)
  if(!Z$no.response || is.null(Z$expr))
    stop("formula should be in the form of x or ~ x | g1*g2*...")
  subset <- eval(substitute(subset), data)
  x <- eval(Z$expr[[1]], data)[subset]
  x <- x[!is.na(x)]
  if(missing(breaks)) {
    if(missing(h)) h <- diff(pretty(x, nbins))[1]
    first <- floor((min(x) - x0)/h)
    last <- ceiling((max(x) - x0)/h)
    breaks <- x0 + h * c(first:last)
  }
  if(any(diff(breaks) <= 0)) stop("breaks must be strictly increasing")
  if(min(x) < min(breaks) || max(x) > max(breaks))
     stop("breaks do not cover the data")
  db <- diff(breaks)
  if(!prob && sqrt(var(db)) > mean(db)/1000)
    warning("Uneven breaks with prob = F will give a misleading plot")
  setup.2d.trellis(formula, data = data, prepanel = prepanel,
		   prepanel.arg = list(breaks = breaks, prob = prob),
		   breaks = breaks, prob = prob, panel = panel,
		   aspect = aspect, xlab = xlab, ylab = ylab,
		   groups = eval(substitute(groups), data), ...,
		   subset = subset)
}
"prepanel.histplot" <-
function(x, y, breaks, prob=TRUE, ...)
{
  data <- x[!is.na(x)]
  x <- breaks
  xl <- max(seq(along=x)[x < min(data)])
  xm <- min(seq(along=x)[x > max(data)])
  x <- x[xl:xm]
  bin <- cut(data, x, include.lowest = TRUE)
  est <- tabulate(bin, length(levels(bin)))
  if (prob) est <- est/(diff(x) * length(data))
  list(xlim = range(x), ylim = c(0, max(est)), dx = diff(x)[-1], dy = diff(est))
}
"panel.histplot" <-
function(x, y, breaks, prob = T, col = trellis.par.get("bar.fill")$col, ...)
{
  data <- x[!is.na(x)]
  x <- breaks
  xl <- max(seq(along=x)[x < min(data)])
  xm <- min(seq(along=x)[x > max(data)])
  x <- x[xl:xm]
  bin <- cut(data, x, include.lowest = TRUE)
  est <- tabulate(bin, length(levels(bin)))
  if (prob) est <- est/(diff(x) * length(data))
  nx <- length(x)
  X <- as.vector(rbind(x[-1], x[ - nx], x[ - nx], x[-1], NA))
  Y <- as.vector(rbind(0, 0, est, est, NA))
  polygon(X, Y, col = col, border = 1, ...)
}



