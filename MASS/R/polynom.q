# file MASS/polynom.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
"polynom"<- function(x, b)
{
    if(!is.loaded(symbol.C("poly"))) dyn.load("horner.o")
    m <- as.integer(length(x))
    n <- as.integer(length(b) - 1)
    storage.mode(x) <- "double"
    p <- x
    storage.mode(b) <- "double"
    .C("poly",
       m,
       val = p,
       x,
       n,
       b)$val
}
