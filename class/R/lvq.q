# file nnet/lvq.q copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
lvqinit <- function(x, cl, size, prior, k=5) {
	x <- as.matrix(x)
	n <- nrow(x)
	p <- ncol(x)
	if(length(cl) != n) stop("x and cl have different lengths")
	g <- as.factor(cl)
	counts <- tapply(rep(1, length(g)), g, sum)
	prop <- counts/n
	np <- length(prop)
# allow for supplied prior
	if(missing(prior)) prior <- prop
	else if(any(prior <0)||round(sum(prior), 5) != 1) 
			stop("invalid prior")
	if(length(prior) != np) stop("prior is of incorrect length")
	if(missing(size)) size <- min(round(0.4 * np * (np-1+p/2),0), n)
	inside <- knn.cv(x, cl, k) == cl
	selected <- numeric(0)
	for(i in 1:np){
		set <- seq(along=g)[unclass(g)==i & inside]
		if(length(set) > 1) 
		   set <- sample(set, min(length(set), round(size*prior[i])))
		selected <- c(selected, set)
	}
	list(x = x[selected, , drop=FALSE], cl = cl[selected])
}

olvq1 <- function(x, cl, codebk, niter = 40*nrow(codebk$x), alpha = 0.3) {
	x <- as.matrix(x)
	n <- nrow(x)
	p <- ncol(x)
	nc <- dim(codebk$x)[1]
	if(length(cl) != n) stop("x and cl have different lengths")
	iters <- sample(n, niter, TRUE)
	z <- .C("VR_olvq",
		as.double(alpha),
		as.integer(n),
		as.integer(p),
		as.double(x),
		as.integer(unclass(cl)),
		as.integer(nc),
		xc = as.double(codebk$x),
		as.integer(codebk$cl),
		as.integer(niter),
		as.integer(iters-1)
		)
	xc <- matrix(z$xc,nc,p)
	dimnames(xc) <- dimnames(codebk$x)
	list(x = xc, cl = codebk$cl)
}

lvq1 <- function(x, cl, codebk, niter = 100*nrow(codebk$x), alpha = 0.03) {
	x <- as.matrix(x)
	n <- nrow(x)
	p <- ncol(x)
	nc <- dim(codebk$x)[1]
	if(length(cl) != n) stop("x and cl have different lengths")
	iters <- sample(n, niter, TRUE)
	z <- .C("VR_lvq1",
		as.double(alpha),
		as.integer(n),
		as.integer(p),
		as.double(x),
		as.integer(unclass(cl)),
		as.integer(nc),
		xc = as.double(codebk$x),
		as.integer(codebk$cl),
		as.integer(niter),
		as.integer(iters-1)
		)
	xc <- matrix(z$xc,nc,p)
	dimnames(xc) <- dimnames(codebk$x)
	list(x = xc, cl = codebk$cl)
}

lvq2 <- function(x, cl, codebk, niter = 100*nrow(codebk$x), alpha = 0.03, win = 0.3) {
	x <- as.matrix(x)
	n <- nrow(x)
	p <- ncol(x)
	nc <- dim(codebk$x)[1]
	if(length(cl) != n) stop("x and cl have different lengths")
	iters <- sample(n, niter, TRUE)
	z <- .C("VR_lvq2",
		as.double(alpha),
		as.double(win),
		as.integer(n),
		as.integer(p),
		as.double(x),
		as.integer(unclass(cl)),
		as.integer(nc),
		xc = as.double(codebk$x),
		as.integer(codebk$cl),
		as.integer(niter),
		as.integer(iters-1)
		)
	xc <- matrix(z$xc,nc,p)
	dimnames(xc) <- dimnames(codebk$x)
	list(x = xc, cl = codebk$cl)
}

lvq3 <- function(x, cl, codebk, niter = 100*nrow(codebk$x), alpha = 0.03, win = 0.3, epsilon = 0.1) {
	x <- as.matrix(x)
	n <- nrow(x)
	p <- ncol(x)
	nc <- dim(codebk$x)[1]
	if(length(cl) != n) stop("x and cl have different lengths")
	iters <- sample(n, niter, TRUE)
	z <- .C("VR_lvq3",
		as.double(alpha),
		as.double(win),
		as.double(epsilon),
		as.integer(n),
		as.integer(p),
		as.double(x),
		as.integer(unclass(cl)),
		as.integer(nc),
		xc = as.double(codebk$x),
		as.integer(codebk$cl),
		as.integer(niter),
		as.integer(iters-1)
		)
	xc <- matrix(z$xc,nc,p)
	dimnames(xc) <- dimnames(codebk$x)
	list(x = xc, cl = codebk$cl)
}

lvqtest <- function(codebk, test) knn1(codebk$x, test, codebk$cl)
