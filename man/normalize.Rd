\name{normalize}
\alias{normalize}
\title{
Global normalization of count data
}
\description{
Normalize a numerical matrix by scaling each column so that the scaled column sums are equal.
}
\usage{
normalize(d)
}
\arguments{
  \item{d}{
A numerical matrix.
}
}
\details{
The average of column sums is computed. A scaling factor is calculated for each colunm so that the scaled column sum is equal to the computed average value.
}
\value{
A matrix of the same size as \code{d} is returned. The column sums of this matrix are equal.
}
\references{
Pham TV (2021). countdata: The Beta-Binomial Test for Count Data. R package version 1.1. https://CRAN.R-project.org/package=countdata
}
\author{
Thang V. Pham
}
\examples{
d <- rbind(c(2.5, 11.2, 7.2),
           c(9.1, 2.2, 7.1))

colSums(d)
# 11.6 13.4 14.3

colSums(normalize(d))
# 13.1 13.1 13.1
}
\keyword{data normalization}
