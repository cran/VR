% file MASS/petrol.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{petrol}
\alias{petrol}
\title{
N. L. Prater's Petrol Refinery Data
}
\section{SUMMARY}{
The yield of a petroleum refining process with four covariates.
The crude oil appears to come from only 10 distinct samples.


These data were originally used by Prater (1956) to
build an estimation equation for the yield of the refining
process of crude oil to gasoline.
}
\section{DATA DESCRIPTION}{
The variables are as follows
}
\arguments{
\item{No}{
Crude oil sample identification label. (factor)
}
\item{SG}{
Specific gravity, degrees API.  (Constant within sample.)
}
\item{VP}{
Vapour pressure in psi. (Constant within sample.)
}
\item{V10}{
Volatility of crude; ASTM 10\% point. (Constant within sample.)
}
\item{EP}{
Desired volatility of gasoline. (The end point.  Varies within sample.)
}
\item{Y}{
Yield as a percentage of crude.
}}
\section{SOURCE}{
N. H. Prater (1956)
Estimate gasoline yields from crudes.
\emph{Petroleum Refiner}
\bold{35}, 236-238.


This dataset is also given in
D. J. Hand, F. Daly, K.  McConway, D.  Lunn, and E. Ostrowski
E. (eds) (1993)
\emph{A Handbook of Small Data Sets.}
Chapman \& Hall.
}
\examples{
library(lme)
data(petrol)
Petrol <- petrol
Petrol[, 2:5] <- scale(as.matrix(Petrol[, 2:5]), scale = F)
pet3.lme <- lme(Y ~ SG + VP + V10 + EP,
                random = ~ 1 | No, data = Petrol, REML=T)  
pet3.lme <- update(pet3.lme, REML=F)
pet4.lme <- update(pet3.lme, fixed = Y ~ V10 + EP)
anova(pet4.lme, pet3.lme)
}
\keyword{datasets}
% Converted by Sd2Rd version 0.3-1.
