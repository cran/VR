% file MASS/Rabbit.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Rabbit}
\alias{Rabbit}
\title{
Blood Pressure in Rabbits
}
\description{
Five rabbits were studied on two occasions, after treatment with
saline (control) and after treatment with the 5-HT_3 antagonist MDL
72222.  After each treatment ascending doses of phenylbiguanide were
injected intravenously at 10 minute intervals and the responses of
mean blood pressure measured.  The goal was to test whether the
cardiogenic chemoreflex elicited by phenylbiguanide depends on the
activation of 5-HT_3 receptors.
}
\format{
This data frame contains 60 rows and the following variables:
}
\arguments{
\item{BPchange}{
change in blood pressure relative to the start of the experiment
}
\item{Dose}{
dose of Phenylbiguanide in micrograms
}
\item{Run}{
label of run (\code{"C1"} to \code{"C5"}, then \code{"M1"} to \code{"M5"})
}
\item{Treatment}{
placebo or the 5-HT_3 antagonist MDL 72222
}
\item{Animal}{
label of animal used (\code{"R1"} to \code{"R5"})
}}
\source{
J. Ludbrook (1994)
Repeated measurements and multiple comparisons in cardiovascular research.
\emph{Cardiovascular Research}
\bold{28}, 303-311.
[The numerical data are not in the paper but were supplied by 
Professor Ludbrook]
}
\keyword{datasets}
