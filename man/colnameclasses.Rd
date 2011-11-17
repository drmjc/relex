\name{colnameclasses}
\alias{colnameclasses}
\title{Print out the colnames and their classes in a 2 column table}
\usage{
  colnameclasses(x)
}
\arguments{
  \item{x}{a table}
}
\value{
  a 2 column table of colname, and colclasses
}
\description{
  Print out the colnames and their classes in a 2 column
  table
}
\examples{
y <- as.data.frame(cbind(1:5, LETTERS[1:5], rep(TRUE, 5)), stringsAsFactors=TRUE)
colnameclasses(y)
}
\author{
  Mark Cowley, 27/3/07
}
\seealso{
  \code{\link{colclasses}}
}

