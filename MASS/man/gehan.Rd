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
\format{
This data frame contains the following columns:
}
\arguments{
\item{pair}{
label for pair
}
\item{time}{
remission time in weeks
}
\item{cens}{
censoring, 0/1
}
\item{treat}{
treatment, control or 6-MP
}}
\source{
Cox, D. R. and Oakes, D. (1984)
\emph{Analysis of Survival Data.}
Chapman \& Hall, p. 7. Taken from


Gehan, E.A. (1965) A generalized Wilcoxon test for comparing 
arbitrarily single-censored samples.
\emph{Biometrika}
\bold{52}, 203-233.
}
\examples{
library(survival4)
data(gehan)
attach(gehan)
gehan.surv <- survfit(Surv(time, cens) ~ treat)
gehan.surv
plot(gehan.surv, conf.int=TRUE, lty=c(3,4), log=TRUE,
   xlab="time of remission (weeks)", ylab="survival")
gehan.exp <- survreg(Surv(time, cens) ~ treat + factor(pair), dist="exp")
summary(gehan.exp)
gehan.exp <- survreg(Surv(time, cens) ~ treat, dist="exp")
summary(gehan.exp)
gehan.wei <- survreg(Surv(time, cens) ~ treat)
summary(gehan.wei)
gehan.cox <- coxph(Surv(time, cens) ~ treat)
summary(gehan.cox)
}
\keyword{datasets}
