% file MASS/petrol.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{petrol}
\alias{petrol}
\title{
N. L. Prater's Petrol Refinery Data
}
\description{
The yield of a petroleum refining process with four covariates.
The crude oil appears to come from only 10 distinct samples.


These data were originally used by Prater (1956) to
build an estimation equation for the yield of the refining
process of crude oil to gasoline.
}
\usage{
data(petrol)
}
\format{
The variables are as follows
\describe{
\item{\code{No}}{
Crude oil sample identification label. (factor)
}
\item{\code{SG}}{
Specific gravity, degrees API.  (Constant within sample.)
}
\item{\code{VP}}{
Vapour pressure in psi. (Constant within sample.)
}
\item{\code{V10}}{
Volatility of crude; ASTM 10\% point. (Constant within sample.)
}
\item{\code{EP}}{
Desired volatility of gasoline. (The end point.  Varies within sample.)
}
\item{\code{Y}}{
Yield as a percentage of crude.
}}}
\source{
N. H. Prater (1956)
Estimate gasoline yields from crudes.
\emph{Petroleum Refiner}
\bold{35}, 236--238.


This dataset is also given in
D. J. Hand, F. Daly, K. McConway, D. Lunn and E. Ostrowski (eds) (1994)
\emph{A Handbook of Small Data Sets.}
Chapman & Hall.
}
\examples{
library(nlme)
Petrol <- petrol
Petrol[, 2:5] <- scale(as.matrix(Petrol[, 2:5]), scale = FALSE)
pet3.lme <- lme(Y ~ SG + VP + V10 + EP,
                random = ~ 1 | No, data = Petrol)
pet3.lme <- update(pet3.lme, method = "ML")
pet4.lme <- update(pet3.lme, fixed = Y ~ V10 + EP)
anova(pet4.lme, pet3.lme)
}
\keyword{datasets}
