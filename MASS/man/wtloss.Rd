% file MASS/wtloss.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{wtloss}
\alias{wtloss}
\title{
Weight Loss Data from an Obese Patient
}
\description{
The data frame gives the weight, in kilograms, of an obese patient at 52
time points over an 8 month period of a weight rehabilitation programme.
}
\usage{
data(wtloss)
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{Days}}{
Time in days since the start of the programme.
}
\item{\code{Weight}}{
Weight in kilograms of the patient.
}}}
\source{
Dr T. Davies, Adelaide.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\examples{
library(nls)
wtloss.fm <- nls(Weight ~ b0 + b1*2^(-Days/th),
    data = wtloss, start = list(b0=90, b1=95, th=120),
    trace = TRUE)
}
\keyword{datasets}
