# file MASS/isoMDS.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
isoMDS <- function(d, y=cmdscale(d, 2), maxit=50, trace=TRUE)
{
    if(any(!is.finite(d))) stop("NAs/Infs not allowed in d")
    if(is.null(n <- attr(d, "Size"))) {
        x <- as.matrix(d)
        if((n <- nrow(x)) != ncol(x))
            stop("Distances must be result of dist or a square matrix")
    } else {
        x <- matrix(0, n, n)
        x[row(x) > col(x)] <- d
        x <- x + t(x)
    }
    if (any(ab <- x[row(x) < col(x)] <= 0)) {
        aa <- cbind(as.vector(row(x)), as.vector(col(x)))[row(x) < col(x),]
        aa <- aa[ab, , drop=FALSE]
        stop(paste("zero or negative distance between objects",aa[1,1],
                   "and", aa[1,2]))
    }
    dis <- x[row(x) > col(x)]
    ord <- order(dis)
    nd <- length(ord)
    n <- nrow(y)
    k <- ncol(y)
    .C("VR_mds_init_data",
       as.integer(nd),
       as.integer(k),
       as.integer(n),
       as.integer(ord - 1),
       as.integer(order(ord) - 1),
       as.double(y)
       )
    tmp <- .C("VR_mds_dovm",
              val = double(1),
              as.integer(maxit),
              as.integer(trace),
              y = as.double(y)
              )
    .C("VR_mds_unload")
    list(points = matrix(tmp$y,,k), stress = tmp$val)
}

Shepard <- function(d, x)
{
#
# Given a dissimilarity d and configuration x, compute Shephard plot
#
  n <- nrow(x)
  k <- ncol(x)
  y <- dist(x)
  ord <- order(d)
  y <- y[ord]
  nd <- length(ord)
  Z <- .C("VR_mds_fn",
	  as.double(y),
	  yf=as.double(y),
	  as.integer(nd),
	  ssq = double(1),
	  as.integer(order(ord)-1),
	  as.double(x),
	  as.integer(n),
	  as.integer(k),
	  g=double(n*k),
	  as.integer(1)
	  )
  list(x = d[ord], y = y, yf = Z$yf)
}

