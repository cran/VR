\name{nlschools}
\alias{nlschools}
\title{
Eighth-Grade Pupils in the Netherlands
}
\description{
Snijders and Bosker (1999) use as a running example a study of 2287
eighth-grade pupils (aged about 11) in 132 classes in 131 schools in
the Netherlands.  Only the variables used in our examples are supplied. 
}
\usage{
data(nlschools)
}
\format{
This data frame contains 2287 rows and the following columns:
\describe{
\item{\code{lang}}{
language test score
}
\item{\code{IQ}}{
Verbal IQ
}
\item{\code{GS}}{
Class size
}
\item{\code{SES}}{
Social-economic status.
}
\item{\code{COMB}}{
were the pupils taught in a multi-grade class (0/1)?
}}}
\source{
Snijders, T. A. B. and Bosker, R. J. (1999) 
\emph{Multilevel Analysis. An Introduction to Basic and Advanced Multilevel Modelling.}
London: Sage.
}
\examples{
data(nlschools)
library(nlme)
nl1 <- nlschools
attach(nl1)
classMeans <- tapply(IQ, class, mean)
nl1$IQave <- classMeans[as.character(class)]
nl1$IQ <- nl1$IQ - nl1$IQave
detach()
cen <- c("IQ", "IQave", "SES")
nl1[cen] <- scale(nl1[cen], center = TRUE, scale = FALSE)

nl.lme <- lme(lang ~ IQ*COMB + IQave + SES, 
              random = ~ IQ | class, data = nl1)
summary(nl.lme)
}
\keyword{datasets}
