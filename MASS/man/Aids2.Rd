% file MASS/Aids2.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Aids2}
\alias{Aids2}
\title{
Australian AIDS Survival Data
}
\description{
Data on patients diagnosed with AIDS in Australia before 1 July 1991.
}
\format{
This data frame contains 2843 rows and the following columns:
}
\arguments{
\item{state}{
Grouped state of origin: NSW includes ACT and \code{"other"} is WA, SA, NT
and TAS.
}
\item{sex}{
Sex of patient
}
\item{diag}{
(Julian) date of diagnosis
}
\item{death}{
(Julian) date of death or end of observation
}
\item{status}{
\code{"A"} (alive) or \code{"D"} (dead) at end of observation
}
\item{T.categ}{
Reported transmission category
}
\item{age}{
Age (years) at diagnosis
}}
\note{
This data set has been slightly jittered as a
condition of its release, to ensure patient confidentiality.
}
\source{
Dr P. J. Solomon and the Australian National Centre in HIV Epidemiology
and Clinical Research.
}
\keyword{datasets}
