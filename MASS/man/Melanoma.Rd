% file MASS/Melanoma.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Melanoma}
\alias{Melanoma}
\title{
Survival from Malignant Melanoma
}
\description{
The \code{Melanoma} data frame has data on 205 patients in Denmark
with malignant melanoma.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{time}{
survival time in days, possibly censored
}
\item{status}{
\code{1} died from melanoma, \code{2} alive, \code{3} dead from other causes
}
\item{sex}{
\code{1} = male, \code{2} = female
}
\item{age}{
age in years
}
\item{year}{
of operation
}
\item{thickness}{
tumour thickness in mm
}
\item{ulcer}{
\code{1} = presence, \code{0} = absence
}}
\source{
P. K. Andersen, O. Borgan, R. D. Gill, and N. Keiding (1993)
\emph{Statistical Models based on Counting Processes.}
Springer.
}
\keyword{datasets}
