# file MASS/stdres.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
lmwork <- function(object)
{
    resid <- object$resid
    hat <- lm.influence(object)$hat
    ok <- !(is.na(resid))
    n.miss <- sum(!ok)
    switch(ifelse(n.miss > 2, 2, n.miss),
           warning("1 missing observation deleted"),
           warning(paste(n.miss, "missing observations deleted")))
    resid <- resid[ok ]
    n <- length(resid)
    p <- object$rank
    rdf <- object$df.resid
    if(is.null(rdf))
        rdf <- n - p
    if(!is.null(object$weights)) {
        wt <- object$weights[ok]
        resid <- resid * wt^0.5
        excl <- wt == 0
        if(any(excl)){
            warning(paste(sum(excl),
                          "rows with zero weights not counted"))
            resid <- resid[!excl]
            if(is.null(object$df.resid))
                rdf <- rdf - sum(excl)
        }
    }
    stdres <- studres <- resid
    if(n > p) {
        stddev <- sqrt(sum(resid^2)/rdf)
        sr <- resid/(sqrt(1 - hat) * stddev)
        stdres <- sr
        studres <- sr/sqrt((n-p-sr^2)/(n-p-1))
    }
    else stddev <- stdres[] <- studres[]<- NA
    list(stdedv=stddev, stdres=stdres, studres=studres)
}
stdres <- function(object) lmwork(object)$stdres
studres <- function(object) lmwork(object)$studres
