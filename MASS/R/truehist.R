# file MASS/truehist.q
# copyright (C) 1994-2002 W. N. Venables and B. D. Ripley
#
"truehist"<-
function(data, nbins = "Scott", h, x0 = -h/1000, breaks,
         prob = TRUE, xlim = range(breaks), ymax = max(est),
	 col = 5,
	 xlab = deparse(substitute(data)), bty = "n", ...)
{
    plot.truehist <-
        function(breaks, est, xlim, ymax, bty, xlab,
                 density = NULL, angle = 45,
                 col = NULL, border = NULL, lty = NULL, lwd = par("lwd"), ...)
    {
        n <- length(breaks)
        plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = "",
             bty = bty, ...)
        rect(breaks[-n], 0, breaks[-1], est,
             density = density, angle = angle,
             col = col, border = border, lty = lty, lwd = lwd)
    }
    xlab  # force evaluation
    data <- data[!is.na(data)]
    if(missing(breaks)) {
        if(missing(h)) {
            if(is.character(nbins))
                nbins <- switch(casefold(nbins),
                                scott = nclass.scott(data),
                                "freedman-diaconis" = , fd = nclass.FD(data)
                                )
            h <- diff(pretty(data, nbins))[1]
        }
        first <- floor((min(data) - x0)/h)
        last <- ceiling((max(data) - x0)/h)
        breaks <- x0 + h * c(first:last)
    }
    if(any(diff(breaks) <= 0)) stop("breaks must be strictly increasing")
    if(min(data) < min(breaks) || max(data) > max(breaks))
        stop("breaks do not cover the data")
    db <- diff(breaks)
    if(!prob && sqrt(var(db)) > mean(db)/1000)
        warning("Uneven breaks with prob = FALSE will give a misleading plot")
    bin <- cut(data, breaks, include.lowest = TRUE)
    est <- tabulate(bin, length(levels(bin)))
    if(prob) est <- est/(diff(breaks) * length(data))
    plot.truehist(breaks, est, xlim, ymax, bty=bty, xlab=xlab,
                  col = col, ...)
#     n <- length(breaks)
#     plot(xlim, c(0, ymax), type = "n", xlab = xlab, ylab = "", bty = bty, ...)
#     rect(breaks[-n], 0, breaks[-1], est, col = col, ...)
    invisible()
}
