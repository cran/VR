% file MASS/rock.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{rock-MASS}
\alias{rock}
\title{
  Measurements on Petroleum Rock Samples
}
\description{
Twelve core samples from petroleum reservoirs were sampled by 4
cross-sections.  Each core sample was measured for permeability, and
each cross-section has total area of pores, total perimeter of pores,
and shape.
}
\usage{
data(rock)
}
\format{
    A data frame with 48 rows and 4 numeric columns.
    \tabular{ll}{
      \code{area}  \tab area of pores space, in pixels
      out of 256 by 256 \cr
     	\code{peri}  \tab perimeter in pixels \cr
      \code{shape} \tab perimeter/sqrt(area) \cr
      \code{perm}  \tab permeability in milli-Darcies
    }
}
\source{
  Data from BP Research, image analysis by Ronit Katz, U. Oxford.

  Katz, R. (1995) \emph{Spatial Analysis of Pore Images}. D.Phil. thesis,
  Department of Statistics, University of Oxford.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\keyword{datasets}
