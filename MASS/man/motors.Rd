% file MASS/motors.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{motorette}
\alias{motorette}
\title{
Accelerated Life Testing of Motorettes
}
\description{
The \code{motorette} data frame has 40 rows and 3 columns.  It describes an
accelerated life test at each of four temperatures of 10 motorettes,
and has rather discrete times.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{temp}{
the temperature (degrees C) of the test
}
\item{time}{
the time in hours to failure or censoring at 8064 hours (= 336 days).
}
\item{cens}{
an indicator variable for death
}}
\source{
Kalbfleisch, J. D. and Prentice, R. L. (1980)
\emph{The Statistical Analysis of Failure Time Data.}
New York: Wiley.


taken from


Nelson, W. D. and Hahn, G. J. (1972)
Linear regression of a regression relationship from censored data.
Part 1 -- simple methods and their application.
\emph{Technometrics }
\bold{14}, 247-276.
}
\examples{
library(survival4)
data(motors)
attach(motors)
plot(survfit(Surv(time, cens) ~ factor(temp)), conf.int=FALSE)
# fit Weibull model
motor.wei <- survreg(Surv(time, cens) ~ temp)
summary(motor.wei)
# and predict at 130C
at130 <- predict(motor.wei, data.frame(temp=130), se.fit=TRUE)
exp(at130$fit)


motor.cox <- coxph(Surv(time, cens) ~ temp)
summary(motor.cox)
# predict at temperature 200
plot( survfit(motor.cox, newdata=data.frame(temp=200)) )
}
\keyword{datasets}
