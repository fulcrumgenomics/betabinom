\name{bb.test}
\alias{bb.test}
\title{
The beta-binomial test
}
\description{
Performs the beta-binomial test for count data.
}
\usage{
bb.test(x, tx, group, alternative = c("two.sided", "less", "greater"),
        n.threads = -1, verbose = TRUE)
}

\arguments{
  \item{x}{
A vector or matrix of counts. When \code{x} a matrix, the test is performed row by row.
}
  \item{tx}{
A vector or matrix of the total sample counts. When \code{tx} is a matrix, the number of rows must be equal to the number of rows of \code{x}.
}
  \item{group}{
A vector of group indicators.
}
  \item{alternative}{
A character string specifying the alternative hypothesis: "two.sided" (default), "greater" or "less".

}
  \item{n.threads}{
The number of threads to be used. When \code{n.threads} is 0, the maximal number of CPU cores is used. When \code{n.threads} is -1 (default), one CPU core less than the maximum is used, and so on.
}
  \item{verbose}{
A logical value. If \code{TRUE} (default), status information is printed.
}

}
\details{
This test is designed for independent samples, two or more groups.
}
\value{
A list with a single component is returned
\item{p.value }{The p-value of the test.}
}
\references{
Pham TV, Piersma SR, Warmoes M, Jimenez CR (2010) On the beta binomial model for analysis of spectral count data in label-free tandem mass spectrometry-based proteomics. \emph{Bioinformatics}, 26(3):363-369.
}
\author{
Thang V. Pham
}
\examples{
x <- c(1, 5, 1, 10, 9, 11, 2, 8)

tx <- c(19609, 19053, 19235, 19374, 18868, 19018, 18844, 19271)

group <- c(rep("cancer", 3), rep("normal", 5))

bb.test(x, tx, group)
# p.value = 0.01568598
}
\keyword{beta-binomial}
\keyword{count data}
\keyword{significance analysis}
\keyword{independent sample test}

