# file MASS/isoMDS.q
# copyright (C) 1994-2001 W. N. Venables and B. D. Ripley
#
isoMDS <- function(d, y = cmdscale(d, k), k = 2, maxit = 50, trace = TRUE)
{
    if(any(!is.finite(as.vector(d)))) stop("NAs/Infs not allowed in d")
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
    if(!is.matrix(y)) stop("y must be a matrix")
    if(any(dim(y) != c(n, k)) ) stop("invalid initial configuration")
    on.exit(.C("VR_mds_unload"))
    .C("VR_mds_init_data",
       as.integer(nd),
       as.integer(k),
       as.integer(n),
       as.integer(ord - 1),
       as.integer(order(ord) - 1),
       as.double(y)
       , PACKAGE = "MASS"
       )
    tmp <- .C("VR_mds_dovm",
              val = double(1),
              as.integer(maxit),
              as.integer(trace),
              y = as.double(y)
              , PACKAGE = "MASS"
              )
    points <- matrix(tmp$y,,k)
    rn <- if(is.matrix(d)) rownames(d) else names(d)
    dimnames(points) <- list(rn, NULL)
    list(points = points, stress = tmp$val)
}

Shepard <- function(d, x)
{
#
# Given a dissimilarity d and configuration x, compute Shepard plot
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
          , PACKAGE = "MASS"
	  )
  list(x = d[ord], y = y, yf = Z$yf)
}

