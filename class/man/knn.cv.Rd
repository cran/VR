% file class/knn.cv.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{knn.cv}
\alias{knn.cv}
\title{
k-Nearest Neighbour Cross-Validatory Classification
}
\description{
k-nearest neighbour cross-validatory classification from training set. For
each row of the training set, the \code{k} nearest (in Euclidean distance) other
training set vectors are found, and the classification is decided by
majority vote, with ties broken at random. If there are ties for the
\code{k}th nearest vector, all candidates are included in the vote.
}
\usage{
knn.cv(train, class, k=1, l=1, prob=FALSE, use.all=TRUE)
}
\arguments{
\item{train}{
matrix or data frame of training set cases.
}
\item{class}{
factor of true classifications of training set
}
\item{k}{
number of neighbours considered.
}
\item{l}{
minimum vote for definite decision, otherwise \code{doubt}. (More
precisely, less than \code{k-l} dissenting votes are allowed, even if \code{k}
is increased by ties.)
}
\item{prob}{
If this is true, the proportion of the votes for the winning class
are returned as attribute \code{prob}.
}
\item{use.all}{
controls handling of ties. If true, all distances equal to the \code{k}th
largest are included. If false, a random selection of distances
equal to the \code{k}th is chosen to use exactly \code{k} neighbours.
}}
\value{
factor of classifications of training set. \code{doubt} will be returned as \code{NA}.
}
\seealso{
\code{\link{knn}}
}
\examples{
data(iris3)
train <- rbind(iris3[1:25,,1],iris3[1:25,,2],iris3[1:25,,3])
test <- rbind(iris3[26:50,,1],iris3[26:50,,2],iris3[26:50,,3])
cl <- factor(c(rep("s",25),rep("c",25), rep("v",25)))
knn.cv(train, cl, k=3, prob=TRUE)
attributes(.Last.value)
}
\keyword{classif}
