% file MASS/caith.d
% copyright (C) 1999 W. N. Venables and B. D. Ripley
%
\name{caith}
\alias{caith}
\title{
Colours of Eyes and Hair of People in Caithness
}
\description{
Data on the cross-classification of people in Caithness, Scotland, by
eye and hair colour. The region of the UK is particularly interesting
as there is a mixture of people of Nordic, Celtic and Anglo-Saxon origin.
}
\usage{
data(caith)
}
\format{
A 4 by 5 table with rows the eye colours (blue, light, medium, dark) and
columns the hair colours (fair, red, medium, dark, black).
}
\source{
Fisher, R.A. (1940) The precision of discriminant functions.
\emph{Annals of Eugenics (London)}
\bold{10}, 422--429.
}
\examples{
corresp(caith)
dimnames(caith)[[2]] <- c("F", "R", "M", "D", "B")
par(mfcol=c(1,3))
library(mva)
plot(corresp(caith, nf=2)); title("symmetric")
plot(corresp(caith, nf=2), type="rows"); title("rows")
plot(corresp(caith, nf=2), type="col"); title("columns")
par(mfrow=c(1,1))
}
\keyword{datasets}
