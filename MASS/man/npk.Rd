% file MASS/npk.d
% copyright (C) 1999 W. N. Venables and B. D. Ripley
%
\name{npk}
\alias{npk}
\title{
Classical N, P, K Factorial Experiment
}
\description{
A classical N, P, K (nitrogen, phosphate, potassium) factorial
experiment on the growth of peas conducted on 6 blocks. Each half of a
fractional factorial design confounding the NPK interaction was used
on 3 of the plots.
}
\format{
The \code{npk} data frame has 24 rows and 5 columns:
}
\arguments{
\item{block}{
which block (label 1 to 6).
}
\item{N}{
indicator (0/1) for the application of nitrogen.
}
\item{P}{
indicator (0/1) for the application of phosphate.
}
\item{K}{
indicator (0/1) for the application of potassium.
}
\item{yield}{
Yield of peas, in pounds/plot (the plots were (1/70) acre).
}}
\source{
Imperial College, London, M.Sc. exercise sheet.
}
\examples{
data(npk)
npk.aov <- aov(yield ~ block + N*P*K, npk)
npk.aov
summary(npk.aov)
alias(npk.aov)
coef(npk.aov)
options(contrasts=c("contr.treatment", "contr.poly"))
npk.aov1 <- aov(yield ~ block + N + K, npk)
summary.lm(npk.aov1)
se.contrast(npk.aov1, list(N=="0", N=="1"), data=npk)
model.tables(npk.aov1, type="means", se=TRUE)
}
\keyword{datasets}
