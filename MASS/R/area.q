# file MASS/area.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
"area"<-
function(f, a, b, ..., fa = f(a, ...), fb = f(b, ...), limit
	 = 10, eps = 1e-5)
{
    h <- b - a
    d <- (a + b)/2
    fd <- f(d, ...)
    a1 <- ((fa + fb) * h)/2
    a2 <- ((fa + 4 * fd + fb) * h)/6
    if(abs(a1 - a2) < eps)
        return(a2)
    if(limit == 0) {
        warning(paste("iteration limit reached near x = ",
                      d))
        return(a2)
    }
    Recall(f, a, d, ..., fa = fa, fb = fd, limit = limit - 1,
           eps = eps) + Recall(f, d, b, ..., fa = fd, fb =
           fb, limit = limit - 1, eps = eps)
}
"fbeta"<-
function(x, alpha, beta)
{
    x^(alpha - 1) * (1 - x)^(beta - 1)
}
"print.abbrev"<-
function(obj, ...)
{
    if(is.list(obj))
        obj <- unlist(obj)
    NextMethod("print")
}
