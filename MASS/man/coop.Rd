% file MASS/coop.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{coop}
\alias{coop}
\title{
Co-operative Trial in Analytical Chemistry
}
\description{
Seven specimens were sent to 6 laboratories in 3 separate batches and
each analysed for Analyte.  Each analysis was duplicated.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{Lab}{
Laboratory, \code{L1}, \code{L2}, \dots{}, \code{L6}.
}
\item{Spc}{
Specimen, \code{S1}, \code{S2}, \dots{}, \code{S6}.
}
\item{Bat}{
Batch, \code{B1}, \code{B2}, \code{B3} (nested within \code{Spc/Lab}),
}
\item{Conc}{
Concentration of Analyte in g/kg.
}}
\source{
Analytical Methods Committee (1987)  
Recommendations for the conduct and
interpretation of co-operative trials, 
\emph{The Analyst}
\bold{112}, 679-686.
}
\seealso{
\code{\link{chem}}, \code{\link{abbey}}.
}
\keyword{datasets}
