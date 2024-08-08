\name{ibb.test}
\alias{ibb.test}
\title{
The inverted beta-binomial test
}
\description{
Performs the inverted beta-binomial test for paired count data.
}
\usage{
ibb.test(x, tx, group, alternative = c("two.sided", "less", "greater"),
         n.threads = -1, BIG = 1e4, verbose = TRUE)
}

\arguments{
  \item{x}{
A vector or matrix of counts. When \code{x} is a matrix, the test is performed row by row.
}
  \item{tx}{
A vector or matrix of the total sample counts. When \code{tx} is a matrix, the number of rows must be equal to the number of rows of \code{x}.
}
  \item{group}{
A vector of group indicators. There should be two groups of equal size. The samples are matched by the order of appearance in each group.
}
  \item{alternative}{
A character string specifying the alternative hypothesis: "two.sided" (default), "greater" or "less".

}
  \item{n.threads}{
The number of threads to be used. When \code{n.threads} is 0, the maximal number of CPU cores is used. When \code{n.threads} is -1 (default), one CPU core less than the maximum is used, and so on.
}
  \item{BIG}{
A number representing a big value of the result, i.e. black-and-white regulation.
}
  \item{verbose}{
A logical value. If \code{TRUE} (default), status information is printed.
}
}
\details{
This test is designed for paired samples, for example data acquired before and after treatment.
}
\value{
A list of values is returned
\item{p.value}{The p-value of the test.}
\item{fc}{An estimation of the common fold change for all sample pairs. A positive value means up-regulation, i.e. the second group is higher, and a negative value down-regulation. A black-and-white regulation is denoted by the \code{BIG} value.}
}
\references{
Pham TV, Jimenez CR (2012) An accurate paired sample test for count data.
\emph{Bioinformatics}, 28(18):i596-i602.
}
\author{
Thang V. Pham
}
\examples{
x <- c(33, 32, 86, 51, 52, 149)

tx <- c(7742608, 15581382, 20933491, 7126839, 13842297, 14760103)

group <- c(rep("cancer", 3), rep("normal", 3))

ibb.test(x, tx, group)
# p.value = 0.004103636
# fc = 2.137632
}
\keyword{inverted beta-binomial test}
\keyword{count data}
\keyword{significance analysis}
\keyword{paired sample test}
