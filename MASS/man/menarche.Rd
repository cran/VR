% file MASS/menarche.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{menarche}
\alias{menarche}
\title{
Age of Menarche data
}
\description{
Proportions of female children at various ages during adolescence
who have reached menarche.
}
\usage{
data(menarche)
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{Age}}{
Average age of the group.  (The groups are reasonably age homogeneous.)
}
\item{\code{Total}}{
Total number of children in the group.
}
\item{\code{Menarche}}{
Number who have reached menarche.
}}}
\source{
Milicer, H. and Szczotka, F. (1966) Age at Menarche in Warsaw girls in
1965.
\emph{Human Biology}
\bold{38}, 199--203.

The data are also given in\cr
Aranda-Ordaz, F.J. (1981)
On two families of transformations to additivity for binary response data.
\emph{Biometrika}
\bold{68}, 357--363.
}
\examples{
data(menarche)
mprob <- glm(cbind(Menarche, Total - Menarche) ~ Age,
             binomial(link = probit), data = menarche)
}
\keyword{datasets}
