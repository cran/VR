# file MASS/write.matrix.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
write.matrix <- function(x, file="", sep=" ")
{
    x <- as.matrix(x)
    p <- ncol(x)
    cat(colnames(x),format(t(x)), file=file, sep=c(rep(sep, p-1), "\n"))
}
