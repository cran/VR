# file MASS/enlist.q
# copyright (C) 1994-8 W. N. Venables and B. D. Ripley
#
"enlist"<-
function(vec)
{
    x <- as.list(vec)
    names(x) <- names(vec)
    x
}
