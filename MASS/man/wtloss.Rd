% file MASS/wtloss.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{wtloss}
\alias{wtloss}
\title{
Weight Loss Data from an Obese Patient
}
\section{SUMMARY}{
The data frame gives the weight, in kilograms, of an obese patient at 52
time points over an 8 month period of a weight rehabilitation programme.
}
\section{DATA DESCRIPTION}{
This data frame contains the following columns:
}
\arguments{
\item{Days}{
Time in days since the start of the programme.
}
\item{Weight}{
Weight in kilograms of the patient.
}}
\section{SOURCE}{
Dr T. Davies, Adelaide.
}
\examples{
### Not usable in R
\dontrun{
wtloss.fm <- nls(Weight ~ b0 + b1*2^(-Days/th),
    data = wtloss, start = list(b0=90, b1=95, th=120),
    trace = TRUE)
}
}
\keyword{datasets}
% Converted by Sd2Rd version 0.3-2.
