% file MASS/janka.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{janka}
\alias{janka}
\title{
Janka Hardness data
}
\description{
Janka hardness is an important structural property of Australian
timbers, but is difficult to measure directly.  It is, however,
directly related to the density of the timber which is relatively
easy to measure.  In an effort to develop a method of predicting Janka
hardness from density, the Janka hardness (here considered the
response) and the density of 36 Australian eucalypt hardwoods was
measured.
}
\format{
The data frame \code{janka} has two components: 
}
\arguments{
\item{Dens}{
the air-dry timber density for the tree, in pound per cubic foot
}
\item{Hard}{
the corresponding Janka hardness, in Janka units.
}}
\details{
A low-degree polynomial regression may be an appropriate model.  A log or
square-root transformation of hardness may also be useful.
}
\source{
E. J. Williams (1959)
\emph{Regression Analysis.}
Wiley, New York. p. 43.
}
\keyword{datasets}
