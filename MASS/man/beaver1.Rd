% file MASS/beaver1.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{beaver1}
\alias{beaver1}
\title{
Body Temperature Series of Beaver 1
}
\description{
The \code{beaver1} data frame has 114 rows and 4 columns on body
temperature measurements at 10 minute intervals.
}
\format{
Reynolds (1994) describes a
small part of a study of the long-term temperature dynamics of beaver
\emph{Castor canadensis}
in north-central Wisconsin.  Body temperature was measured by
telemetry every 10 minutes for four females, but data from a one
period of less than a day for each of two animals is used there.


This data frame contains the following columns:
}
\arguments{
\item{day}{
Day of observation (in days since the beginning of 1990),
December 12-13.
}
\item{time}{
Time of observation, in the form \code{0330} for 3.30am
}
\item{temp}{
Measured body temperature in degrees Celcius
}
\item{activ}{
Indicator of activity outside the retreat
}}
\section{NOTE}{
The observation at 22:20 is missing.
}
\source{
P. S. Reynolds (1994) Time-series analyses of beaver body temperatures.
Chapter 11 of
Lange, N., Ryan, L., Billard, L., Brillinger, D., Conquest, L. 
and Greenhouse, J. eds (1994)
\emph{Case Studies in Biometry.}
New York: John Wiley and Sons.
}
\seealso{
\code{\link{beaver2}}
}
\examples{
\dotnrun{### Not usable in R
attach(beav1)
beav1$hours <- 24*(day-346) + trunc(time/100) + (time\%\%100)/60
detach()
plot(beav1$hours, beav1$temp, type="l", xlab="time", 
   ylab="temperature", main="Beaver 1")
usr <- par("usr"); usr[3:4] <- c(-0.2, 8); par(usr=usr)
lines(beav1$hours, beav1$activ, type="s", lty=2)
temp <- rts(c(temp[1:82], NA, temp[83:114]), start=9.5, 
            frequency=6, units="hours")
activ <- rts(c(activ[1:82], NA, activ[83:114]), start=9.5,
             frequency=6, units="hours")


acf(temp[1:53]) # and also type="partial"
ar(temp[1:53])
act <- c(rep(0, 10), activ)
X <- cbind(1, act=act[11:125], act1 = act[10:124], 
          act2 = act[9:123], act3 = act[8:122])
arima.mle(temp, xreg=X, model=list(ar=0.82))
}}
\keyword{datasets}
