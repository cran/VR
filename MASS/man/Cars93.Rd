% file MASS/Cars93.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Cars93}
\alias{Cars93}
\title{
Data from 93 Cars on Sale in the USA in 1993
}
\description{
The \code{Cars93} data frame has 93 rows and 27 columns.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{Manufacturer}{
Manufacturer
}
\item{Model}{
Model
}
\item{Type}{
Type: Small, Sporty, Compact, Midsize, Large, Van
}
\item{Min.Price}{
Minimum Price (in \$1,000) - price for a basic version
}
\item{Price}{
Midrange Price (in \$1,000) - average of \code{Min.Price} and \code{Max.Price}
}
\item{Max.Price}{
Maximum Price (in \$1,000) - price for "a premium version"
}
\item{MPG.city}{
City MPG (miles per US gallon by EPA rating)
}
\item{MPG.highway}{
Highway MPG
}
\item{AirBags}{
Air Bags standard. Factor: none, driver only, or driver & passenger
}
\item{DriveTrain}{
Drive train type: rear wheel, front wheel or 4WD; (factor).
}
\item{Cylinders}{
Number of cylinders (missing for Mazda RX-7, which has a rotary engine).
}
\item{EngineSize}{
Engine size (litres)
}
\item{Horsepower}{
Horsepower (maximum)
}
\item{RPM}{
RPM (revs per minute at maximum horsepower)
}
\item{Rev.per.mile}{
Engine revolutions per mile (in highest gear)
}
\item{Man.trans.avail}{
Is a manual transmission version available? (yes or no, Factor).
}
\item{Fuel.tank.capacity}{
Fuel tank capacity (US gallons)
}
\item{Passengers}{
Passenger capacity (persons)
}
\item{Length}{
Length  (inches)
}
\item{Wheelbase}{
Wheelbase (inches)
}
\item{Width}{
Width (inches)
}
\item{Turn.circle}{
U-turn space (feet)
}
\item{Rear.seat.room}{
Rear seat room (inches) (missing for 2-seater vehicles)
}
\item{Luggage.room}{
Luggage capacity (cubic feet) (missing for vans)
}
\item{Weight}{
Weight (pounds)
}
\item{Origin}{
Of non-USA or USA company origins? (factor)
}
\item{Make}{
Combination of Manufacturer and Model (character)
}}
\source{
Cars were selected at random from among 1993 passenger car models that
were listed in both the \emph{Consumer Reports} issue and the 
\emph{PACE Buying Guide}.  Pickup trucks and Sport/Utility vehicles were
eliminated due to incomplete information in the \emph{Consumer Reports}
source.  Duplicate models (e.g., Dodge Shadow and Plymouth Sundance)
were listed at most once.


Further description can be found in Lock (1993). Use the URL
\url{http://jse.stat.ncsu.edu/jse/v1n1/datasets.lock}
or send the message
     \url{send jse/v1n1/datasets.lock}
to the email address \url{archive@jse.stat.ncsu.edu}
}
\source{
Lock, R. H. (1993)
1993 New Car Data.
\emph{Journal of Statistics Education}
\bold{1}(1)
}
\keyword{datasets}
