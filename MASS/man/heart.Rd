\name{heart}
\alias{heart}
\title{
  The Stanford heart transplant data
}
\description{
  The Stanford heart transplant data from Kalbfliesch & Prentice.
  A data frame.
}
\details{
The components are:
}
\arguments{
\item{start}{
start time in days. Either 0 or time of transplant.
}
\item{transplant}{
is this before (0) or after (1) transplantation?
}
\item{stop}{
time of transplant or death or censoring (days)
}
\item{event}{
death (1) or censored (0)
}
\item{age}{
patient's age in years at acceptance, less 48.
}
\item{year}{
date of acceptance. Actually days since 1967/10/1 divided by 365.25.
}
\item{surgery}{
  prior surgery (0 = no, 1 = yes).
}
\item{id}{
  Patient's ID
}}
\source{
Kalbfleisch, J.D. and Prentice R.L. (1980)
\emph{The Statistical Analysis of Failure Time Data.} 
}
\keyword{datasets}
