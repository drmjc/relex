\name{unionN}
\alias{unionN}
\title{Union of multiple vectors...}
\usage{unionN(...)
}
\description{Union of multiple vectors
This function extends the \code{\link{union}} function to allow more than 2
vectors to be unioned. In practice, we do not use union at all, but a more
optimised algorithm.}
\value{A vector containing the union of all the input vectors.}
\author{Mark Cowley}
\seealso{\code{\link{intersectN}}, \code{\link{union}}}
\keyword{logic}
\keyword{array}
\arguments{\item{\dots}{at least 2 vectors to be unioned}
}
\examples{a <- letters[1:15]
b <- letters[5:20]
c <- letters[10:25]
unionN(a, b, c)}

