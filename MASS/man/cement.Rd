% file MASS/cement.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{cement}
\alias{cement}
\title{
Heat Evolved by Setting Cements
}
\description{
Experiment on the heat evolved in the setting of each of 13
cements.
}
\format{
}
\arguments{
\item{"x1, x2, x3, x4"}{
Proportions (\%) of active ingredients
}
\item{y}{
heat evolved in cals/gm
}}
\source{
Hald, A. (1957)
\emph{Statistical Theory with Engineering Applications.}
Wiley, New York.
}
\examples{
data(cement)
lm(y ~ x1 + x2 + x3 + x4, cement)
}
\keyword{datasets}
