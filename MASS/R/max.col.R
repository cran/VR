# file MASS/max.col.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
max.col <- function(m)
{
    if(!is.loaded(symbol.C("VR_max_col")))
        stop("Compiled code has not been dynamically loaded")
    m <- as.matrix(m)
    n <- nrow(m)
    .C("VR_max_col",
       as.double(m),
       as.integer(n),
       as.integer(ncol(m)),
       rmax = integer(n),
       NAOK=TRUE
       )$rmax
}

