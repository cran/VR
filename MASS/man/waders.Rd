% file MASS/waders.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{waders}
\alias{waders}
\title{
Counts of Waders at 15 Sites in South Africa.
}
\description{
The \code{waders} data frame has 15 rows and 19 columns.
The entries are counts of waders in summer.
}
\format{
This data frame contains the following columns (species)
}
\arguments{
\item{S1}{
Oystercatcher
}
\item{S2}{
White-fronted Plover
}
\item{S3}{
Kitt Lutz's Plover
}
\item{S4}{
Three-banded Plover
}
\item{S5}{
Grey Plover
}
\item{S6}{
Ringed Plover
}
\item{S7}{
Bar-tailed Godwit
}
\item{S8}{
Whimbrel
}
\item{S9}{
Marsh Sandpiper
}
\item{S10}{
Greenshank
}
\item{S11}{
Common Sandpiper
}
\item{S12}{
Turnstone
}
\item{S13}{
Knot
}
\item{S14}{
Sanderling
}
\item{S15}{
Little Stint
}
\item{S16}{
Curlew Sandpiper
}
\item{S17}{
Ruff
}
\item{S18}{
Avocet
}
\item{S19}{
Black-winged Stilt


The rows are the sites:
A = Namibia North coast
B = Namibia North wetland
C = Namibia South coast
D = Namibia South wetland
E = Cape North coast
F = Cape North wetland
G = Cape West coast
H = Cape West wetland
I = Cape South coast
J= Cape South wetland
K = Cape East coast
L = Cape East wetland
M = Transkei coast
N = Natal coast
O = Natal wetland
}}
\source{
J.C. Gower and D.J. Hand (1996)
\emph{Biplots}
Chapman \& Hall Table 9.1. Quoted as from:


R.W. Summers, L.G. Underhill, D.J. Pearson and D.A. Scott (1987)
Wader migration systems in south and eastern Africa and western Asia.
\emph{Wader Study Group Bulletin}
\bold{49} Supplement, 15-34.
}
\examples{
library(mva)
data(waders)
plot(corresp(waders, nf=2))
}
\keyword{datasets}
