% file MASS/whiteside.d
% copyright (C) 1999 W. N. Venables and B. D. Ripley
%
\name{whiteside}
\alias{whiteside}
\title{
House Insulation: Whiteside's Data
}
\description{
Mr Derek Whiteside of the UK Building Research Station recorded the
weekly gas consumption and average external temperature at his own
house in south-east England for 26 weeks before, and 30 weeks after
cavity-wall insulation was installed. The object of the exercise was
to assess the effect of the insulation on gas consumption.
}
\format{
The \code{whiteside} data frame has 56 rows and 3 columns.:
}
\arguments{
\item{Insul}{
A factor, before or after insulation.
}
\item{Temp}{
Purportedly the average outside temperature in degrees Celsius. (These
values is far too low for any 56-week period in the 1960s in
South-East England. It might be the weekly average of daily minima.)
}
\item{Gas}{
The weekly gas consumption in 1000s of cubic feet.
}}
\source{
A data set collected in the 1960s by Mr Derek Whiteside of the
UK~Building Research Station. Reported by


Hand, D. J., Daly, F., McConway, K., Lunn, D. and Ostrowski, E. eds (1993)
\emph{A Handbook of Small Data Sets.}
Chapman \& Hall, p. 69.
}
\examples{
data(whiteside)
coplot(Gas ~ Temp | Insul, whiteside, panel = 
  function(x, y, \dots) {
     points(x, y, \dots)
     abline(lm(y ~ x))
  }, xlab = "Average external temperature (deg. C)",
  ylab = "Gas consumption  (1000 cubic feet)")


gasB <- lm(Gas ~ Temp, whiteside, subset=Insul=="Before")
gasA <- update(gasB, subset=Insul=="After")
summary(gasB)
summary(gasA)
gasBA <- lm(Gas ~ Insul/Temp - 1, whiteside)
summary(gasBA)


gasQ <- lm(Gas ~ Insul/(Temp + I(Temp^2)) - 1, whiteside)
summary(gasQ)$coef


gasPR <- lm(Gas ~ Insul + Temp, whiteside)  
anova(gasPR, gasBA)
options(contrasts = c("contr.treatment", "contr.poly"))
gasBA1 <- lm(Gas ~ Insul*Temp, whiteside)
summary(gasBA1)$coef
}
\keyword{datasets}
