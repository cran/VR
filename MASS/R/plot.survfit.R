# based on code in survival4 and modifications in MASS library.
plot.survfit <-
function(surv, conf.int, mark.time = TRUE, mark = 3, col = 1, lty = 1,
         lwd = 1, cex = 1, logplot = FALSE, yscale = 1, xscale = 1,
         xlab = "", ylab = "", xaxs = "i", add = FALSE, cloglog=FALSE,
         strata, ...)
{
    put.line <- function(xx, yy, cloglog, lty, col, lwd)
    {
        nn <- length(xx)
        if (nn <= 1) return()
        if (cloglog) yy <- -log(yy)[-1]
        whom <- c(match(unique(yy[-nn]), yy), nn)
        lines(xx[whom], yy[whom], type="s", lty=lty, col=col, lwd=lwd)
    }

    if(!inherits(surv, "survfit"))
        stop("First arg must be the result of survfit")

    stime <- surv$time / xscale
    ssurv <- surv$surv
    if(missing(conf.int)) {
        if(is.null(surv$strata) && !is.matrix(ssurv)) conf.int <- TRUE
        else conf.int <- FALSE
    }
    if(is.null(surv$strata)) {
        nstrat <- 1
        stemp <- rep(1, length(surv$time))
    }
    else {
        nstrat <- length(surv$strata)
        stemp <- rep(1:nstrat, surv$strata)
    }
    if(is.null(surv$n.event)) mark.time <- FALSE
    # expected survival curve
    # set default values for missing parameters
    if (is.matrix(ssurv)) ncurve <- nstrat * ncol(ssurv)
    else                  ncurve <- nstrat
    mark <- rep(mark, length = ncurve)
    col <- rep(col, length = ncurve)
    lty <- rep(lty, length = ncurve)
    lwd <- rep(lwd, length = ncurve)
    if(!is.logical(mark.time) && is.numeric(mark.time))
        mark.time <- sort(mark.time[mark.time > 0])
    if (missing(xaxs)) temp <- 1.04*max(stime)
    else               temp <- max(stime)

    #
    # for log plots we have to be tricky about the y axis scaling
    #
    if(!add) {
        if(cloglog) {
            yy <- ssurv[ssurv > 0 & ssurv < 1]
            # strategy is to force 10-20% on-scale
            ymin <- min(0.1, yy, na.rm=TRUE)
            ymax <- max(0.2, yy, na.rm=TRUE)
            plot(c(min(stime), temp),
                 yscale * c(-log(ymax), -log(ymin)), type = "n",
                 log = "xy", xlab = xlab, ylab = ylab, xaxs = xaxs, ...)
        }
        else if(logplot) {
            ymin <- min(0.1, ssurv[!is.na(ssurv) & ssurv > 0])
            ssurv[!is.na(ssurv) & ssurv == 0] <- ymin
            plot(c(0, temp), yscale * c(0.99, ymin), type = "n",
                 log = "y", xlab = xlab, ylab = ylab, xaxs = xaxs, ...)
        }
        else plot(c(0, temp), yscale * c(0, 1), type = "n",
                  xlab = xlab, ylab = ylab, xaxs = xaxs, ...)
    }
    if(yscale != 1) par(usr = par("usr")/c(1, 1, yscale, yscale))
    #
    # put up the curves one by one
    #
    i <- 0
    xend <- NULL
    yend <- NULL
    if(missing(strata)) strata <- unique(stemp)
    for(j in strata) {
	who <- (stemp == j)
	if (cloglog) xx <- stime[who] else xx <- c(0, stime[who])
	nn <- length(xx)
	deaths <- c(-1, surv$n.event[who])
	if (is.matrix(ssurv)) {
	    for (k in 1:ncol(ssurv)) {
		i <- i + 1
		yy <- c(1, ssurv[who,k])
		ind <- deaths == 0
                if (cloglog) {
                    yy <- -log(yy)[-1]
                    ind <- ind[-1]
                }
		put.line(xx, yy, FALSE, lty[i], col[i], lwd[i])
		if (is.numeric(mark.time)) {
		    indx <- mark.time
		    for (k in seq(along=mark.time))
			indx[k] <- sum(mark.time[k] > xx)
		    points(mark.time[indx < nn], yy[indx[indx < nn]],
			   pch=mark[i], col=col[i], cex=cex)
		} else if (mark.time == TRUE && any(ind))
		    points(xx[ind], yy[ind], pch=mark[i], col=col[i], cex=cex)
		xend <- c(xend,max(xx))
		yend <- c(yend,min(yy))

		if (conf.int == TRUE && !is.null(surv$upper)) {
		    if (ncurve == 1) lty[i] <- lty[i] + 1
		    yy <- c(1,surv$upper[who,k])
		    put.line(xx, yy, cloglog, lty[i], col[i], lwd[i])
		    yy <- c(1,surv$lower[who,k])
		    put.line(xx, yy, cloglog, lty[i], col[i], lwd[i])
		}
	    }
	} else {
            i <- i + 1
            yy <- c(1, ssurv[who])
            ind <- deaths == 0
            if (cloglog) {
                yy <- -log(yy)[-1]
                ind <- ind[-1]
            }
            put.line(xx, yy, FALSE, lty[i], col[i], lwd[i])
            if(is.numeric(mark.time)) {
                nn <- length(xx)
                indx <- mark.time
                for(k in seq(along = mark.time))
                    indx[k] <- sum(mark.time[k] > xx)
                points(mark.time[indx < nn], yy[indx[indx < nn]], pch
                       = mark[i], col = col[i], cex = cex)
            }
            else if(mark.time == TRUE && any(ind))
                points(xx[ind], yy[ind], pch = mark[i], col = col[i], cex = cex)
            xend <- c(xend, max(xx))
            yend <- c(yend, min(yy))
            if(conf.int == TRUE && !is.null(surv$upper)) {
                if(ncurve == 1) lty[i] <- lty[i] + 1
                yy <- c(1, surv$upper[who])
                put.line(xx, yy, cloglog, lty[i], col[i], lwd[i])
                yy <- c(1, surv$lower[who])
                put.line(xx, yy, cloglog, lty[i], col[i], lwd[i])
            }
        }
    }
    invisible(list(x = xend, y = yend))
}
