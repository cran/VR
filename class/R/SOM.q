batchSOM <- function(data, grid = somgrid(), radii, init)
{
    data <- as.matrix(data)
    nd <- nrow(data)
    ng <- nrow(grid$pts)
    if(missing(init))
        init <- data[sample(1:nd, ng, replace = FALSE),]
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
    eqscplot(c(0, x$xdim+(x$topo == "hexagonal") + 1),
             c(x$ydim + 1, 0),
             axes = F, type = "n", xlab = "", ylab = "", ...)
    if(type == "p") points(x$pts, cex = 2, ...)
    invisible()
}

plot.SOM <- function(x, ...)
{
    if(!inherits(x, "SOM")) stop("wrong plot method used")
    eqscplot(c(0, x$grid$xdim+(x$grid$topo == "hexagonal") + 1),
             c(x$grid$ydim + 1, 0),
             axes = F, type = "n", xlab = "", ylab = "", ...)
    stars(x$codes, location = x$grid$pts, labels = NULL, len = 0.5)
    invisible()
}
