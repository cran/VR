% file MASS/gehan.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{gehan}
\alias{gehan}
\title{
Remission Times of Leukaemia Patients
}
\description{
A data frame from a trial of 42 leukaemia patients. Some were
treated with the drug 
\emph{6-mercaptopurine}
and the rest are controls.  The trial was designed as matched pairs,
both withdrawn from the trial when either came out of remission.
}
\usage{
data(gehan)
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{pair}}{
label for pair
}
\item{\code{time}}{
remission time in weeks
}
\item{\code{cens}}{
censoring, 0/1
}
\item{\code{treat}}{
treatment, control or 6-MP
}}}
\source{
Cox, D. R. and Oakes, D. (1984)
\emph{Analysis of Survival Data.}
Chapman & Hall, p. 7. Taken from


Gehan, E.A. (1965) A generalized Wilcoxon test for comparing 
arbitrarily single-censored samples.
\emph{Biometrika}
\bold{52}, 203--233.
}
\examples{
library(survival)
data(gehan)
gehan.surv <- survfit(Surv(time, cens) ~ treat, data=gehan,
     conf.type="log-log")
summary(gehan.surv)
survreg(Surv(time, cens) ~ factor(pair) + treat, gehan, dist="exp")
summary(survreg(Surv(time, cens) ~ treat, gehan, dist="exp"))
summary(survreg(Surv(time, cens) ~ treat, gehan))
gehan.cox <- coxph(Surv(time, cens) ~ treat, gehan)
summary(gehan.cox)
}
\keyword{datasets}
