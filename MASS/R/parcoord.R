parcoord <- function(x, col = 1, lty = 1, ...)
{
    x <- apply(x, 2, function(x) (x - min(x))/(max(x) - min(x)))
    matplot(1:ncol(x), t(x), type = "l", col = col, lty = lty,
            xlab="", ylab = "",
            axes = FALSE, ...)
    axis(1, at = 1:ncol(x), labels = colnames(x))
    for(i in 1:ncol(x)) lines(c(i, i), c(0, 1), col = "grey70")
    invisible()
}
