batchSOM <- function(data, grid = somgrid(), radii, init)
{
    data <- as.matrix(data)
    nd <- nrow(data)
    ng <- nrow(grid$pts)
    if(missing(init))
        init <- data[sample(1:nd, ng, replace = FALSE), , drop = FALSE]
    nhbrdist <- as.matrix(dist(grid$pts))
    for(r in radii) {
        cl <- as.numeric(knn1(init, data, 1:ng))
        A <- (nhbrdist <= r)[, cl]
        ind <- rowSums(A) > 0
        init[ind, ] <- A[ind, ] %*% data / rowSums(A)[ind]
    }
    structure(list(grid = grid, codes = init), class = "SOM")
}

somgrid <- function(xdim = 8, ydim = 6, topo = c("rectangular", "hexagonal"))
{
    topo <- match.arg(topo)
    x <- 1:xdim
    y <- 1:ydim
    pts <- as.matrix(expand.grid(x = x, y = y))
    if(topo == "hexagonal") {
       pts[, 1] <- pts[, 1] + 0.5 * (pts[, 2] %% 2)
       pts[, 2] <- sqrt(3)/2 * pts[, 2]
    }
    res <- list(pts = pts, xdim = xdim, ydim = ydim, topo = topo)
    class(res) <- "somgrid"
    res
}

plot.somgrid <- function(x, type = "p", ...)
{
    if(!inherits(x, "somgrid")) stop("wrong plot method used")
    MASS::eqscplot(c(0, x$xdim+(x$topo == "hexagonal") + 1),
                   c(x$ydim + 1, 0),
                   axes = FALSE, type = "n", xlab = "", ylab = "", ...)
    if(type == "p") points(x$pts, cex = 2, ...)
    invisible()
}

plot.SOM <- function(x, ...)
{
    if(!inherits(x, "SOM")) stop("wrong plot method used")
    MASS::eqscplot(c(0, x$grid$xdim+(x$grid$topo == "hexagonal") + 1),
                   c(x$grid$ydim + 1, 0),
                   axes = FALSE, type = "n", xlab = "", ylab = "", ...)
    stars(x$codes, location = x$grid$pts, labels = NULL, len = 0.5)
    invisible()
}

SOM <- function(data, grid = somgrid(), rlen = 10000,
                alpha = seq(0.05, 0, len = rlen),
                radii = seq(4, 1, len = rlen),
                init)
{
    data <- as.matrix(data)
    nd <- nrow(data)
    if(!nd) stop("SOM called with no data")
    ng <- nrow(grid$pts)
    nphases <- 1
    if(is.list(alpha)) {
        nphases <- length(alpha)
        if(!is.list(radii) || length(radii) != nphases)
            stop("'radii' must be a list of the same length as 'alpha'")
    }
    if(missing(init))
        init <- data[sample(1:nd, ng, replace = FALSE), , drop = FALSE]
    codes <- init
    nhbrdist <- as.matrix(dist(grid$pts))
    if(nphases == 1) {
        rlen <- length(alpha)
        if(length(radii) != rlen) stop("alpha and radii do not match")
        codes <- .C("VR_onlineSOM",
                    data = as.double(data),
                    codes = as.double(codes),
                    nhbrdist = as.double(nhbrdist),
                    alpha = as.double(alpha),
                    radii = as.double(radii),
                    n = as.integer(nrow(data)),
                    p = as.integer(ncol(data)),
                    ncodes = as.integer(nrow(init)),
                    rlen = as.integer(rlen)
                    , PACKAGE = "class"
                    )$codes
    } else {
        for(k in 1:nphases) {
            rlen <- length(alpha[[k]])
            if(length(radii[[k]]) != rlen)
                stop("alpha and radii do not match")
            codes <- .C("VR_onlineSOM",
                        data = as.double(data),
                        codes = as.double(codes),
                        nhbrdist = as.double(nhbrdist),
                        alpha = as.double(alpha[[k]]),
                        radii = as.double(radii[[k]]),
                        n = as.integer(nrow(data)),
                        p = as.integer(ncol(data)),
                        ncodes = as.integer(nrow(init)),
                        rlen = as.integer(rlen)
                        , PACKAGE = "class"
                        )$codes
        }
    }
    dim(codes) <- dim(init)
    colnames(codes) <- colnames(init)
    structure(list(grid = grid, codes = codes), class = "SOM")
}
