% file MASS/birthwt.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{birthwt}
\alias{birthwt}
\title{
Risk Factors Associated with Low Infant Birth Weight
}
\description{
The \code{birthwt} data frame has 189 rows and 10 columns.
The data were collected at Baystate Medical Center, Springfield, Mass
during 1986.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{low}{
indicator of birth weight less than 2.5kg
}
\item{age}{
mother's age in years
}
\item{lwt}{
mother's weight in pounds at last menstrual period
}
\item{race}{
mother's race (\code{1} = white, \code{2} = black, \code{3} = other)
}
\item{smoke}{
smoking status during pregnancy
}
\item{ptl}{
number of previous premature labours
}
\item{ht}{
history of hypertension
}
\item{ui}{
presence of uterine irritability
}
\item{ftv}{
number of physician visits during the first trimester
}
\item{bwt}{
birth weight in grams
}}
\source{
Hosmer, D.W. and Lemeshow, S. (1989)
\emph{Applied Logistic Regression.}
New York: Wiley
}
\examples{
data(birthwt)
attach(birthwt)
race <- factor(race, labels=c("white", "black", "other"))
ptd <- factor(ptl > 0)
ftv <- factor(ftv)
levels(ftv)[-(1:2)] <- "2+"
bwt <- data.frame(low=factor(low), age, lwt, race,
    smoke=(smoke>0), ptd, ht=(ht>0), ui=(ui>0), ftv)
detach("birthwt")
options(contrasts = c("contr.treatment", "contr.poly"))
glm(low ~ ., binomial, bwt)
}
\keyword{datasets}
