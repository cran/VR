% file MASS/max.col.d copyright (C) 1999 W. N. Venables and B. D. Ripley
%
\name{max.col}
\alias{max.col}
\title{
Find Maximum Position in Matrix
}
\description{
Find the maximum position for each row of a matrix, breaking ties at random.
}
\usage{
max.col(A)
}
\arguments{
\item{A}{
a matrix
}}
\value{
index of a maximal value for each row.
}
\details{
Ties are broken at random. The determination of 'tie' assumes that the
entries are probabilities.
}
\keyword{utilities}
