# file nnet/knn.q copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
knn1 <- function(train, test, cl)
{
	train <- as.matrix(train)
	if(is.null(dim(test))) dim(test) <- c(1, length(test))
	test <- as.matrix(test)
	p <- ncol(train)
	ntr <- nrow(train)
	if(length(cl) != ntr) stop("train and class have different lengths")
	nte <- nrow(test)
	if(ncol(test) != p) stop("Dims of test and train differ")
	clf <- as.factor(cl)
	nc <- max(unclass(clf))
	res <- .C("VR_knn1",
		as.integer(ntr),
		as.integer(nte),
		as.integer(p),
		as.double(train),
		as.integer(unclass(clf)),
		as.double(test),
		res = integer(nte),
		integer(nc+1),
		as.integer(nc),
		d = double(nte)
		)$res
	factor(res, levels=seq(along=levels(clf)), labels=levels(clf))
}

knn <- function(train, test, cl, k=1, l=0, prob=FALSE, use.all=TRUE)
{
	train <- as.matrix(train)
	if(is.null(dim(test))) dim(test) <- c(1, length(test))
	test <- as.matrix(test)
	p <- ncol(train)
	ntr <- nrow(train)
	if(length(cl) != ntr) stop("train and class have different lengths")
	if(ntr < k) {
	   warning(paste("k =",k,"exceeds number",ntr,"of patterns"))
	   k <- ntr
	}
	if (k < 1) stop(paste("k =",k,"must be at least 1"))
	nte <- nrow(test)
	if(ncol(test) != p) stop("Dims of test and train differ")
	clf <- as.factor(cl)
	nc <- max(unclass(clf))
	Z <- .C("VR_knn",
		as.integer(k),
		as.integer(l),
		as.integer(ntr),
		as.integer(nte),
		as.integer(p),
		as.double(train),
		as.integer(unclass(clf)),
		as.double(test),
		res = integer(nte),
		pr = double(nte),
		integer(nc+1),
		as.integer(nc),
		as.integer(FALSE),
		as.integer(use.all)
		)
	res <- factor(Z$res, levels=seq(along=levels(clf)),labels=levels(clf))
	if(prob) attr(res, "prob") <- Z$pr
	res
}

knn.cv <- function(train, cl, k=1, l=0, prob=FALSE, use.all=TRUE)
{
	train <- as.matrix(train)
	p <- ncol(train)
	ntr <- nrow(train)
	if(ntr-1 < k) {
	   warning(paste("k =",k,"exceeds number",ntr-1,"of patterns"))
	   k <- ntr - 1
	}
	if (k < 1) stop(paste("k =",k,"must be at least 1"))
	clf <- as.factor(cl)
	nc <- max(unclass(clf))
	Z <- .C("VR_knn",
		as.integer(k),
		as.integer(l),
		as.integer(ntr),
		as.integer(ntr),
		as.integer(p),
		as.double(train),
		as.integer(unclass(clf)),
		as.double(train),
		res = integer(ntr),
		pr = double(ntr),
		integer(nc+1),
		as.integer(nc),
		as.integer(TRUE),
		as.integer(use.all)
		)
	res <- factor(Z$res, levels=seq(along=levels(clf)),labels=levels(clf))
	if(prob) attr(res, "prob") <- Z$pr
	res
}
