% file MASS/cats.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{cats}
\alias{cats}
\title{
Anatomical Data from Domestic Cats
}
\description{
The heart and body weights of samples of male and female cats used for
digitalis experiments.  The cats were all adult, over 2 kg body weight.
}
\usage{
data(cats)
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{Sex}}{
Sex factor. Levels \code{"Female"} and \code{"Male"}.
}
\item{\code{Bwt}}{
Body weight in kg.
}
\item{\code{Hwt}}{
Heart weight in g.
}}}
\source{
R. A. Fisher (1947) The analysis of covariance method for the relation
between a part and the whole, 
\emph{Biometrics}
\bold{3}, 65--68.
}
\keyword{datasets}
