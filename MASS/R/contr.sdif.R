# file MASS/contr.sdif.q
# copyright (C) 1994-9 W. N. Venables and B. D. Ripley
#
contr.sdif <- function(n, contrasts = TRUE)
{
    # contrasts generator giving 'successive difference' contrasts.
    if(is.numeric(n) && length(n) == 1) {
        if(n %% 1 || n < 2)
            stop("invalid degree")
        lab <- as.character(seq(n))
    } else {
        lab <- as.character(n)
        n <- length(n)
        if(n < 2)
            stop("invalid number of levels")
    }
    if(contrasts) {
        contr <- col(matrix(nrow = n, ncol = n - 1))
        upper.tri <- !lower.tri(contr)
        contr[upper.tri] <- contr[upper.tri] - n
        structure(contr/n,
                  dimnames = list(lab, paste(lab[-1], lab[-n], sep="-")))
    } else structure(diag(n), dimnames = list(lab, lab))
}
