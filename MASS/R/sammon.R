# file MASS/sammon.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
sammon <- function(d, y, k=2, niter=100, trace=TRUE, magic=0.2, tol=1e-4)
{
   call <- match.call()
   if(any(!is.finite(as.vector(d))))
      stop("NAs/Infs not allowed in d")
   if(is.null(n <- attr(d, "Size"))) {
      x <- as.matrix(d)
      if((n <- nrow(x)) != ncol(x))
        stop("Distances must be result of dist or a square matrix")
   }
   else {
      x <- matrix(0, n, n)
      x[row(x) > col(x)] <- d
      x <- x + t(x)
   }
   if (any(ab <- x[row(x) < col(x)]<=0)) {
	aa <- cbind(as.vector(row(x)), as.vector(col(x)))[row(x) < col(x),]
	aa <- aa[ab,,drop=FALSE]
	stop(paste("zero or negative distance between objects", aa[1,1],
	 "and", aa[1,2]))
	}
   if(missing(y)) y <- cmdscale(d, k=k)
   if(any(dim(y) != c(n, k)) ) stop("invalid initial configuration")
   storage.mode(x) <- "double"
   storage.mode(y) <- "double"
   if(!is.loaded(symbol.C("VR_sammon")))
     stop("Compiled code has not been dynamically loaded")
   z <- .C("VR_sammon",
      x = x,
      as.integer(n),
      as.integer(k),
      y = y,
      as.integer(niter),
      e = double(1),
      as.integer(trace),
      as.double(magic),
      as.double(tol)
      )
   list(points=z$y, stress=z$e, call=call)
}

