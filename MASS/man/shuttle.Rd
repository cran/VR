% file MASS/shuttle.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{shuttle}
\alias{shuttle}
\title{
Space Shuttle Autolander Problem
}
\description{
The \code{shuttle} data frame has 256 rows and 7 columns.
The first six columns are categorical variables giving example
conditions; the seventh is the decision. The first 253 rows are the
training set, the last 3 the test conditions.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{stability}{
Stable positioning or not (\code{stab / xstab})
}
\item{error}{
Size of error (\code{MM / SS / LX / XL})
}
\item{sign}{
Sign of error, positive or negative (\code{pp / nn})
}
\item{wind}{
Wind sign (\code{head / tail})
}
\item{magn}{
Wind strength (\code{Light / Medium / Strong / Out of Range})
}
\item{vis}{
Visibility (\code{yes / no})
}
\item{use}{
Use the autolander or not
}}
\source{
D. Michie (1989)
Problems of computer-aided concept formation. In 
\emph{Applications of Expert Systems 2,}
ed. J. R. Quinlan, Turing Institute Press / Addison-Wesley, pp. 310-333.
}
\keyword{datasets}
