\name{heart}
\alias{heart}
\title{
  The Stanford heart transplant data
}
\description{
  The Stanford heart transplant data from Kalbfleisch & Prentice.
}
\format{
  A data frame with columns:
\describe{
\item{\code{start}}{
start time in days. Either 0 or time of transplant.
}
\item{\code{transplant}}{
is this before (0) or after (1) transplantation?
}
\item{\code{stop}}{
time of transplant or death or censoring (days)
}
\item{\code{event}}{
death (1) or censored (0)
}
\item{\code{age}}{
patient's age in years at acceptance, less 48.
}
\item{\code{year}}{
date of acceptance. Actually days since 1967/10/1 divided by 365.25.
}
\item{\code{surgery}}{
  prior surgery (0 = no, 1 = yes).
}
\item{id}{
  Patient's ID
}}
}
\source{
Kalbfleisch, J.D. and Prentice R.L. (1980)
\emph{The Statistical Analysis of Failure Time Data.} 
}
\keyword{datasets}
