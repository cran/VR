# file nnet/nnet.q copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
nnet <- function(object, ...)
{
  if(is.null(class(object))) class(object) <- data.class(object)
  UseMethod("nnet")
}

nnet.formula <- function(formula, data = NULL, weights, ...,
			subset, na.action = na.fail, contrasts=NULL)
{
  class.ind <- function(cl)
  {
    n <- length(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (as.vector(unclass(cl)) - 1)] <- 1
    dimnames(x) <- list(names(cl), levels(cl))
    x
  }
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, sys.frame(sys.parent()))))
      m$data <- as.data.frame(data)
  m$... <- m$contrasts <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.frame(sys.parent()))
  Terms <- attr(m, "terms")
  x <- model.matrix(Terms, m, contrasts)
  xint <- match("(Intercept)", colnames(x), nomatch=0)
  if(xint > 0) x <- x[, -xint, drop=FALSE] # Bias term is used for intercepts
  w <- model.extract(m, weights)
  if(length(w) == 0) w <- rep(1, nrow(x))
  y <- model.extract(m, response)
  if(is.factor(y)) {
    lev <- levels(y)
    counts <- table(y)
    if(any(counts == 0)) {
      warning(paste("group(s)", paste(lev[counts == 0], collapse=" "),
		    "are empty"))
      y <- factor(y, levels=lev[counts > 0])
    }
    if(length(lev) == 2) {
      y <- as.vector(unclass(y)) - 1
      res <- nnet.default(x, y, w, entropy=TRUE, ...)
      res$lev <- lev
    } else {
      y <- class.ind(y)
      res <- nnet.default(x, y, w, softmax=TRUE, ...)
      res$lev <- lev
    }
  } else res <- nnet.default(x, y, w, ...)
  res$terms <- Terms
  res$coefnames <- colnames(x)
  res$call <- match.call()
  if(!is.null(attr(m, "na.action"))) res$na.action <- attr(m, "na.action")
  class(res) <- c("nnet.formula", "nnet")
  res
}

nnet.default <-
function(x, y, weights, size, Wts, mask=rep(TRUE, length(wts)),
	 linout=FALSE, entropy=FALSE, softmax=FALSE, censored=FALSE, skip=FALSE,
	 rang=0.7, decay=0, maxit=100, Hess=FALSE, trace=TRUE,
         MaxNWts=1000, abstol=1.0e-4, reltol=1.0e-8, ...)
{
  net <- NULL
  x <- as.matrix(x)
  y <- as.matrix(y)
  if(any(is.na(x))) stop("missing values in x")
  if(any(is.na(y))) stop("missing values in y")
  if(dim(x)[1] != dim(y)[1]) stop("nrows of x and y must match")
  if(linout && entropy) stop("entropy fit only for logistic units")
  if(softmax) {
    linout <- TRUE
    entropy <- FALSE
  }
  if(censored) {
    linout <- TRUE
    entropy <- FALSE
    softmax <- TRUE
  }
  net$n <- c(dim(x)[2], size, dim(y)[2])
  net$nunits <- 1 + sum(net$n)
  net$nconn <- rep(0, net$nunits+1)
  net$conn <- numeric(0)
  net <- norm.net(net)
  if(skip) net <- add.net(net, seq(1,net$n[1]),
			  seq(1+net$n[1]+net$n[2], net$nunits-1))
  if((nwts <- length(net$conn))==0) stop("No weights to fit")
  if(nwts > MaxNWts)
    stop(paste("Too many (", nwts, ") weights", sep=""))
  nsunits <- net$nunits
  if(linout) nsunits <- net$nunits - net$n[3]
  net$nsunits <- nsunits
  net$decay <- decay
  net$entropy <- entropy
  net$softmax <- softmax
  net$censored <- censored
  if(missing(Wts))
    if(rang > 0) wts <- runif(nwts, -rang, rang)
    else wts <- rep(0, nwts)
  else wts <- Wts
  if(length(wts) != nwts) stop("weights vector of incorrect length")
  if(length(mask) != length(wts)) stop("incorrect length of mask")
  if(trace) {
    cat("# weights: ", length(wts))
    nw <- sum(mask != 0)
    if(nw < length(wts)) cat(" (", nw, " variable)\n",sep="")
    else cat("\n")
  }
  if(length(decay) == 1) decay <- rep(decay, length(wts))
  if(!is.loaded(symbol.C("VR_set_net")))
    stop("Compiled code has not been dynamically loaded")
  .C("VR_set_net",
     as.integer(net$n),
     as.integer(net$nconn),
     as.integer(net$conn),
     as.double(decay),
     as.integer(nsunits),
     as.integer(entropy),
     as.integer(softmax),
     as.integer(censored)
     )
  ntr <- dim(x)[1]
  nout <- dim(y)[2]
  if(missing(weights)) weights <- rep(1, ntr)
  if(length(weights) != ntr || any(weights < 0))
    stop("invalid weights vector")
  Z <- as.double(cbind(x,y))
  storage.mode(weights) <- "double"
  z <- .C("VR_set_train", as.integer(ntr), Z, weights)
  on.exit(.C("VR_unset_train"))
  tmp <- .C("VR_dovm",
	    as.integer(length(wts)),
	    wts=as.double(wts),
	    val=double(1),
	    as.integer(maxit),
	    as.logical(trace),
	    as.integer(mask),
            as.double(abstol), as.double(reltol)
	    )
  net$value <- tmp$val
  net$wts <- tmp$wts
  tmp <- matrix(.C("VR_nntest",
		   as.integer(ntr), Z, tclass = double(ntr*nout),
		   as.double(net$wts))$tclass,  ntr, nout)
  dimnames(tmp) <- list(rownames(x), colnames(y))
  net$fitted.values <- tmp
  tmp <- y - tmp
  dimnames(tmp) <- list(rownames(x), colnames(y))
  net$residuals <- tmp
  .C("VR_unset_net")
  if(entropy) net$lev <- c("0","1")
  if(softmax) net$lev <- colnames(y)
  net$call <- match.call()
  if(Hess) net$Hessian <- nnet.Hess(net, x, y, weights)
  class(net) <- "nnet"
  net
}


predict.nnet <- function(object, newdata, type=c("raw","class"), ...)
{
  if(!inherits(object, "nnet")) stop("object not of class nnet")
  type <- match.arg(type)
  if(missing(newdata)) z <- fitted(object)
  else {
    if(inherits(object, "nnet.formula")) {#
      # formula fit
      newdata <- as.data.frame(newdata)
      rn <- row.names(newdata)
# work hard to predict NA for rows with missing data
      Terms <- delete.response(object$terms)
      m <- model.frame(Terms, newdata, na.action = na.omit)
      keep <- match(row.names(m), rn)
      x <- model.matrix(Terms, m)
      xint <- match("(Intercept)", colnames(x), nomatch=0)
      if(xint > 0) x <- x[, -xint, drop=FALSE] # Bias term is used for intercepts
    } else { #
      # matrix ...  fit
      if(is.null(dim(newdata)))
        dim(newdata) <- c(1, length(newdata))# a row vector
      x <- as.matrix(newdata)		# to cope with dataframes
      if(any(is.na(x))) stop("missing values in x")
      keep <- 1:nrow(x)
      rn <- rownames(x)
    }
    ntr <- nrow(x)
    nout <- object$n[3]
    if(!is.loaded(symbol.C("VR_set_net")))
      stop("Compiled code has not been dynamically loaded")
    .C("VR_set_net",
       as.integer(object$n), as.integer(object$nconn),
       as.integer(object$conn), as.double(object$decay),
       as.integer(object$nsunits), as.integer(0),
       as.integer(object$softmax), as.integer(object$censored))
    z <- matrix(NA, length(keep), nout,
                dimnames = list(rn, dimnames(object$fitted)[[2]]))
    z[keep, ] <- matrix(.C("VR_nntest",
                           as.integer(ntr),
                           as.double(x),
                           tclass = double(ntr*nout),
                           as.double(object$wts)
                        )$tclass, ntr, nout)
    .C("VR_unset_net")
  }
  switch(type, raw = z,
	 class = {
	   if(is.null(object$lev)) stop("inappropriate fit for class")
	   if(ncol(z) > 1) object$lev[max.col(z)]
	   else object$lev[1 + (z > 0.5)]
	 })
}

eval.nn <- function(wts)
{
  z <- .C("VR_dfunc",
	 as.double(wts), df = double(length(wts)), fp = as.double(1))
  fp <- z$fp
  attr(fp, "gradient") <- z$df
  fp
}

add.net <- function(net, from, to)
{
  nconn <- net$nconn
  conn <- net$conn
  for(i in to){
    ns <- nconn[i+2]
    cadd <- from
    if(nconn[i+1] == ns) cadd <- c(0,from)
    con <- NULL
    if(ns > 1) con <- conn[1:ns]
    con <- c(con, cadd)
    if(length(conn) > ns) con <- c(con, conn[(ns+1):length(conn)])
    for(j in (i+1):net$nunits) nconn[j+1] <- nconn[j+1]+length(cadd)
    conn <- con
  }
  net$nconn <- nconn
  net$conn <- con
  net
}

norm.net <- function(net)
{
  n <- net$n; n0 <- n[1]; n1 <- n0+n[2]; n2 <- n1+n[3];
  if(n[2] <= 0) return(net)
  net <- add.net(net, 1:n0,(n0+1):n1)
  add.net(net, (n0+1):n1, (n1+1):n2)
}

which.is.max <- function(x)
{
  y <- seq(along=x)[x == max(x)]
  if(length(y) > 1) sample(y,1) else y
}

nnet.Hess <- function(net, x, y, weights)
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  if(dim(x)[1] != dim(y)[1]) stop("dims of x and y must match")
  nw <- length(net$wts)
  decay <- net$decay
  if(length(decay) == 1) decay <- rep(decay, nw)
  if(!is.loaded(symbol.C("VR_set_net")))
    stop("Compiled code has not been dynamically loaded")
  .C("VR_set_net",
     as.integer(net$n),
     as.integer(net$nconn),
     as.integer(net$conn),
     as.double(decay),
     as.integer(net$nsunits),
     as.integer(net$entropy),
     as.integer(net$softmax),
     as.integer(net$censored)
     )
  ntr <- dim(x)[1]
  nout <- dim(y)[2]
  if(missing(weights)) weights <- rep(1, ntr)
  if(length(weights) != ntr || any(weights < 0))
    stop("invalid weights vector")
  Z <- as.double(cbind(x,y))
  storage.mode(weights) <- "double"
  z <- .C("VR_set_train",
	  as.integer(ntr),
	  Z,
	  weights
	  )
  z <- matrix(.C("VR_nnHessian", as.double(net$wts), H = double(nw*nw))$H,
              nw, nw)
  .C("VR_unset_train"); .C("VR_unset_net")
  z
}

class.ind <- function(cl)
{
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)) )
  x[(1:n) + n*(unclass(cl)-1)] <- 1
  dimnames(x) <- list(names(cl), levels(cl))
  x
}

print.nnet <- function(x, ...)
{
  if(!inherits(x, "nnet")) stop("Not legitimate a neural net fit")
  cat("a ",x$n[1],"-",x$n[2],"-",x$n[3]," network", sep="")
  cat(" with", length(x$wts),"weights\n")
  if(length(x$coefnames))  cat("inputs:", x$coefnames, "\noutput(s):",
                               deparse(formula(x)[[2]]), "\n")
  cat("options were -")
  tconn <- diff(x$nconn)
  if(tconn[length(tconn)] > x$n[2]+1) cat(" skip-layer connections ")
  if(x$nunits > x$nsunits && !x$softmax) cat(" linear output units ")
  if(x$entropy) cat(" entropy fitting ")
  if(x$softmax) cat(" softmax modelling ")
  if(x$decay[1] > 0) cat(" decay=", x$decay[1], sep="")
  cat("\n")
  invisible(x)
}

coef.nnet <- function(object, ...)
{
  wts <- object$wts
  wm <- c("b", paste("i", seq(length=object$n[1]), sep=""))
  if(object$n[2] > 0)
  wm <- c(wm, paste("h", seq(length=object$n[2]), sep=""))
  if(object$n[3] > 1)  wm <- c(wm,
	  paste("o", seq(length=object$n[3]), sep=""))
  else wm <- c(wm, "o")
  names(wts) <- apply(cbind(wm[1+object$conn],
                            wm[1+rep(1:object$nunits - 1, diff(object$nconn))]),
                      1, function(x)  paste(x, collapse = "->"))
  wts
}

summary.nnet <- function(object, ...)
{
  class(object) <- c("summary.nnet", class(object))
  object
}


print.summary.nnet <- function(x, ...)
{
  cat("a ",x$n[1],"-",x$n[2],"-",x$n[3]," network", sep="")
  cat(" with", length(x$wts),"weights\n")
  cat("options were -")
  tconn <- diff(x$nconn)
  if(tconn[length(tconn)] > x$n[2]+1) cat(" skip-layer connections ")
  if(x$nunits > x$nsunits && !x$softmax) cat(" linear output units ")
  if(x$entropy) cat(" entropy fitting ")
  if(x$softmax) cat(" softmax modelling ")
  if(x$decay[1] > 0) cat(" decay=", x$decay[1], sep="")
  cat("\n")
  wts <- format(round(coef.nnet(x),2))
  lapply(split(wts,rep(1:x$nunits, tconn)),
	 function(x) print(x, quote=FALSE))
  invisible(x)
}

residuals.nnet <- function(object, ...) object$residuals
