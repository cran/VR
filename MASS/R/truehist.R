# file MASS/truehist.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
"truehist"<-
function(data, nbins = nclass.scott(data), h, x0 = -h/1000, breaks, prob = TRUE,
	 xlim = range(breaks), ymax = max(est),
	 col = 5,
	 xlab = deparse(substitute(data)), bty = "n", ...)
{
  eval(xlab)
  data <- data[!is.na(data)]
  if(missing(breaks)) {
    if(missing(h)) h <- diff(pretty(data, nbins))[1]
    first <- floor((min(data) - x0)/h)
    last <- ceiling((max(data) - x0)/h)
    breaks <- x0 + h * c(first:last)
  }
  if(any(diff(breaks) <= 0)) stop("breaks must be strictly increasing")
  if(min(data) < min(breaks) || max(data) > max(breaks))
     stop("breaks do not cover the data")
  db <- diff(breaks)
  if(!prob && sqrt(var(db)) > mean(db)/1000)
    warning("Uneven breaks with prob = F will give a misleading plot")
  bin <- cut(data, breaks, include.lowest = TRUE)
  est <- tabulate(bin, length(levels(bin)))
  if(prob) est <- est/(diff(breaks) * length(data))
  n <- length(breaks)
  plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = "", bty = bty)
  rect(breaks[-n], 0, breaks[-1], est, col = col, ...)
  invisible()
}
