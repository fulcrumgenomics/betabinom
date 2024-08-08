\name{fold.change}
\alias{fold.change}
\title{
Fold change calculation
}
\description{
Calculates the fold changes between two numerical matrices row by row.
}
\usage{
fold.change(d1, d2, BIG = 1e4)
}

\arguments{
  \item{d1}{
The first data matrix.
}
  \item{d2}{
The second data matrix.
}
  \item{BIG}{
A number representing a big value of the result, i.e. black-and-white regulation.
}
}
\details{
The two matrices \code{d1} and \code{d2} must have the same number of rows.
A positive value means up-regulation where the average of \code{d2} is higher than that of \code{d1}. Conversely, a negative value means down-regulation where the the average of \code{d1} is higher than that of \code{d2}. If one group contains all zeros, a positive or negative \code{BIG} value is returned.
}
\value{
A vector of fold changes is returned.
}
\references{
Pham TV (2021). countdata: The Beta-Binomial Test for Count Data. R package version 1.1. https://CRAN.R-project.org/package=countdata
}
\author{
Thang V. Pham
}
\examples{
x <- rbind(c(1.5, 1.2, 10.2),
           c(9.1, 2.2, 8.1))

fold.change(x, 2 * x)
# returns a column vector of 2

fold.change(x, 0.5 * x)
# returns a column vector of -2
}
\keyword{fold change}
