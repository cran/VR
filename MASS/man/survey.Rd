% file MASS/survey.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{survey}
\alias{survey}
\title{
Student Survey Data
}
\description{
This data frame contains the responses of 237 Statistics I students at
the University of Adelaide to a number of questions.
}
\usage{
survey
}
\format{
The components of the data frame are:
\describe{
\item{\code{Sex}}{
The sex of the student. (Factor with levels \code{"Male"} and \code{"Female"}.)
}
\item{\code{Wr.Hnd}}{
The span (distance from tip of thumb to tip of little finger of spread
hand) of writing hand, in centimetres.
}
\item{\code{NW.Hnd}}{
Span of non-writing hand
}
\item{\code{W.Hnd}}{
Writing hand of student. (Factor, with levels \code{"Left"} and \code{"Right"}.)
}
\item{\code{Fold}}{
  \dQuote{Fold your arms! Which is on top} (Factor, with levels
  \code{"R on L"}, \code{"L on R"}, \code{"Neither"}.)
}
\item{\code{Pulse}}{
Pulse rate of student (beats per minute.)
}
\item{\code{Clap}}{
  \sQuote{Clap your hands! Which hand is on top?} (Factor, with levels
  \code{"Right"}, \code{"Left"}, \code{"Neither"}.)
}
\item{\code{Exer}}{
How often the student exercises. (Factor, with levels \code{"Freq"}
(frequently), \code{"Some"}, "\code{None"})
}
\item{\code{Smoke}}{
How much the student smokes. (Factor, levels \code{"Heavy"}, \code{"Regul"}
(regularly), \code{"Occas"} (occasionally), \code{"Never"}.)
}
\item{\code{Height}}{
The height of the student, in centimetres.
}
\item{\code{M.I}}{
Indicates whether the student expressed height in imperial
(feet/inches) or metric (centimetres/metres) units. (Factor, levels
\code{"Metric"}, \code{"Imperial"}.)
}
\item{\code{Age}}{
Age of the student, in years.
}}}
\references{
  Venables, W. N. and Ripley, B. D. (1999)
  \emph{Modern Applied Statistics with S-PLUS.} Third
  Edition. Springer.
}
\keyword{datasets}
