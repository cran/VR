# file MASS/corresp.q
# copyright (C) 1994-2004 W. N. Venables and B. D. Ripley
#
corresp <- function(x, ...) UseMethod("corresp")

corresp.xtabs <- function(x, ...)
{
  if((m <- length(dim(x))) > 2)
    stop("frequency table is ", m, "-dimensional")
  corresp.matrix(x, ...)
}

corresp.data.frame <- function(x, ...)
    corresp.matrix(as.matrix(x), ...)

corresp.default <- function(x, ...)
    stop("invalid table specification")

corresp.factor <- function(x, y, ...)
    corresp.matrix(table(x, y), ...)

corresp.formula <- function(formula, data = parent.frame(), ...)
{
    rhs <- formula[[length(formula)]]
    if(length(rhs[[2]]) > 1 || length(rhs[[3]]) > 1)
        stop("higher-way table requested.  Only 2-way allowed")
    tab <- table(eval(rhs[[2]], data), eval(rhs[[3]], data))
    names(dimnames(tab)) <- as.character(c(rhs[[2]], rhs[[3]]))
    corresp.matrix(tab, ...)
}

corresp.matrix <- function(x, nf = 1, ...)
{
    if(any(x < 0 | x %% 1 > 10 * sqrt(.Machine$double.eps)))
        warning("negative or non-integer entries in table")
    if((N <- sum(x)) == 0) stop("all frequencies are zero")
    Dr <- drop(x %*% (rep(1/N, ncol(x))))
    Dc <- drop((rep(1/N, nrow(x))) %*% x)
    if(any(Dr == 0) || any(Dc == 0)) stop("empty row or column in table")
    x1 <- x/N - outer(Dr, Dc)
    Dr <- 1/sqrt(Dr)
    Dc <- 1/sqrt(Dc)
    if(is.null(dimnames(x)))
        dimnames(x) <- list(Row = paste("R", 1:nrow(x)),
                            Col = paste("C", 1:ncol(x)))
    if(is.null(names(dimnames(x))))
        names(dimnames(x)) <- c("Row", "Column")
    X.svd <- svd(t(t(x1 * Dr) * Dc))
    dimnames(X.svd$u) <- list(rownames(x), NULL)
    dimnames(X.svd$v) <- list(colnames(x), NULL)
    res <- list(cor = X.svd$d[1:nf], rscore = X.svd$u[, 1:nf] * Dr,
                cscore = X.svd$v[, 1:nf] * Dc, Freq = x)
    class(res) <- "correspondence"
    res
}

plot.correspondence <- function(x, scale=1, ...)
{
    if(length(x$cor) > 1) return(invisible(biplot(x, ...)))
    Fr <- x$Freq
    rs <- x$rscore
    cs <- x$cscore
    xs <- range(cs)
    xs <- xs + diff(xs) * c(-1/5, 1/5)
    ys <- range(rs)
    ys <- ys + diff(ys) * c(-1/5, 1/5)
    x <- cs[col(Fr)]
    y <- rs[row(Fr)]
    rcn <- names(dimnames(Fr))
    plot(x, y, xlim = xs, ylim = ys, xlab = rcn[2], ylab = rcn[1], pch = 3)
    size <- min(par("pin"))/20 * scale
    symbols(x, y, circles = as.vector(sqrt(Fr)), inches = size, add = TRUE)
    x0 <- (min(cs) + min(xs))/2
    y0 <- (min(rs) + min(ys))/2
    text(cs, y0, names(cs))
    text(x0, rs, names(rs), adj = 1)
    invisible()
}

print.correspondence <- function(x, ...)
{
    cat("First canonical correlation(s):", format(x$cor, ...), "\n")
    rcn <- names(dimnames(x$Freq))
    cat("\n", rcn[1], "scores:\n")
    print(x$rscore)
    cat("\n", rcn[2], "scores:\n")
    print(x$cscore)
    invisible(x)
}

biplot.correspondence <-
    function(x, type = c("symmetric", "rows", "columns"), ...)
{
    if(length(x$cor) < 2) stop("biplot is only possible if nf >= 2")
    type <- match.arg(type)
    X <- x$rscore[, 1:2]
    if(type != "columns") X <- X %*% diag(x$cor[1:2])
    colnames(X) <- rep("", 2)
    Y <- x$cscore[, 1:2]
    if(type != "rows")  Y <- Y %*% diag(x$cor[1:2])
    colnames(Y) <- rep("", 2)
    switch(type, "symmetric" = biplot(X, Y, var.axes = FALSE, ...),
           "rows" = biplot.bdr(X, Y, ...),
           "columns" = biplot.bdr(Y, X, ...))
    points(0, 0, pch = 3, cex = 3)
    invisible()
}

biplot.bdr <-
    function(obs, bivars, col, cex = rep(par("cex"), 2),
             olab = NULL, vlab = NULL, xlim = NULL, ylim = NULL, ...)
{
  # for cases where we need equal scales for the two sets of vars.
    expand.range <- function(x)
    {
        if(x[1] > 0) x[1] <-  - x[1]
        else if(x[2] < 0) x[2] <-  - x[2]
        x
    }
    n <- dim(obs)[1]
    p <- dim(bivars)[1]
    vlab.real <- rownames(bivars)
    if(is.logical(vlab)) vlab <- vlab.real[vlab]
    else if(length(vlab) != p) vlab <- vlab.real
    else vlab <- as.character(vlab)
    if(!length(vlab)) {
        vlab.real <- vlab <- paste("Var", 1:p)
        dimnames(bivars) <- list(vlab, colnames(bivars))
    }
    if(length(olab)) olab <- rep(as.character(olab), length.out = n)
    else {
        olab <- rownames(obs)
        if(length(olab) != n) olab <- as.character(1:n)
    }
    if(length(cex) != 2) cex <- rep(cex, length.out = 2)
    if(missing(col)) {
        col <- par("col")
        if (!is.numeric(col)) col <- match(col, palette())
        col <- c(col, col + 1)
    }
    else if(length(col) != 2) col <- rep(col, length.out = 2)
    ro1 <- expand.range(range(obs[, 1]))
    ro2 <- expand.range(range(obs[, 2]))
    rv1 <- expand.range(range(bivars[, 1]))
    rv2 <- expand.range(range(bivars[, 2]))
    if(!(length(xlim) || length(ylim)))
        xlim <- ylim <- range(ro1, ro2, rv1, rv2)
    else if(!length(xlim)) xlim <- range(ro1, rv1)
    else if(!length(ylim)) ylim <- range(ro2, rv2)
    on.exit(par(oldpar))
    oldpar <- par(pty = "s")
    plot(obs, type = "n", xlim = xlim, ylim = ylim, col = col[1], ...)
    text(obs, labels=olab, cex = cex[1], col = col[1], ...)
    par(new = TRUE)
    plot(bivars, axes = FALSE, type = "n", xlim = xlim, ylim =
         ylim, xlab = "", ylab = "", col = col[1], ...)
    axis(3, col = col[2])
    axis(4, col = col[2])
    box(col = col[1])
    text(bivars, labels=vlab, cex = cex[2], col = col[2], ...) #
    invisible()
}

