# file MASS/write.matrix.q
# copyright (C) 1994-2001 W. N. Venables and B. D. Ripley
#
write.matrix <- function(x, file = "", sep = " ", blocksize)
{
    x <- as.matrix(x)
    p <- ncol(x)
    cn <- colnames(x)
    if(!missing(blocksize) && blocksize > 0) {
        cat(cn, file=file, sep=c(rep(sep, p-1), "\n"))
        nlines <- 0
        nr <- nrow(x)
        while (nlines < nr) {
            nb <- min(blocksize, nr - nlines)
            cat(format(t(x[nlines + (1:nb), ])),
                file = file, append = TRUE,
                sep = c(rep(sep, p-1), "\n"))
            nlines <- nlines + nb
        }
    } else
        cat(c(cn, format(t(x))), file=file, sep=c(rep(sep, p-1), "\n"))
}
