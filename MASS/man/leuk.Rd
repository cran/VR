% file MASS/leuk.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{leuk}
\alias{leuk}
\title{
Survival Times and White Blood Counts for Leukaemia Patients
}
\description{
A data frame of data from 33 leukaemia patients.
}
\usage{
data(leuk)
}
\format{
A data frame with columns:


\describe{
\item{\code{wbc}}{
white blood count
}
\item{\code{ag}}{
a test result, \code{"present"} or \code{"absent"}
}
\item{\code{time}}{
survival time in weeks
}}}
\source{
  Cox, D. R. and Oakes, D. (1984) \emph{Analysis of Survival Data}.
  Chapman & Hall, p. 9.

Taken from

Feigl, P. & Zelen, M. (1965) Estimation of exponential survival
probabilities with concomitant information.
\emph{Biometrics} \bold{21}, 826--838.
}
\examples{
library(survival)
data(leuk)
plot(survfit(Surv(time) ~ ag, data = leuk), lty = 2:3, col = 2:3)

# now Cox models
leuk.cox <- coxph(Surv(time) ~ ag + log(wbc), leuk)
summary(leuk.cox)
}
\keyword{datasets}
