% file MASS/farms.d
% copyright (C) 1998-9 W. N. Venables and B. D. Ripley
%
\name{farms}
\alias{farms}
\title{
Ecological Factors in Farm Management
}
\description{
The \code{farms} data frame has 20 rows and 4 columns. The rows are farms
on the Dutch island of Terschelling and the columns are factors
describing the management of grassland. 
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{Mois}{
Five levels of soil moisture -- level 3 does not occur at these 20 farms.
}
\item{Manag}{
Grassland management type (\code{SF} = standard, \code{BF} = biological, \code{HF} =
hobby farming, \code{NM} = nature conservation).
}
\item{Use}{
Grassland use (\code{U1} = hay production, \code{U2} = intermediate, \code{U3} = grazing)
}
\item{Manure}{
Manure usage -- classes C0 to C4.
}}
\source{
J.C. Gower and D.J. Hand (1996)
\emph{Biplots}
Chapman & Hall, Table 4.6. 


Quoted as from:
R.H.G. Jongman, C.J.F. ter Braak and O.F.R. van Tongeren (1987)
\emph{Data Analysis in Community and Landscape Ecology.}
PUDOC, Wageningen
}
\examples{
data(farms)
farms.mca <- mca(farms, abbrev=T)  # Use levels as names
eqscplot(farms.mca$cs, type="n")
text(farms.mca$rs, cex=0.7)
text(farms.mca$cs, labels=dimnames(farms.mca$cs)[[1]], cex=0.7)
}
\keyword{datasets}
