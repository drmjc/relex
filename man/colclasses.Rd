\name{colclasses}
\alias{colclasses}
\title{Get the column classes of matrix-like data
For matrix-like data, this returns the classes of each column. It's quite
useful for checking whether factors snuck into a \code{data.frame}.  All
columns in a \code{matrix} have the same type, but \code{data.frame}'s may
have mixed column types.}
\usage{
  colclasses(x)
}
\arguments{
  \item{x}{a \code{matrix}, or a \code{data.frame}}
}
\value{
  colclasses: a character vector of column classes, with
  the same length as the number of columns. Typical classes
  include \code{numeric}, \code{character}, \code{logical},
  \code{factor}, etc\dots{}
}
\description{
  Get the column classes of matrix-like data For
  matrix-like data, this returns the classes of each
  column. It's quite useful for checking whether factors
  snuck into a \code{data.frame}.  All columns in a
  \code{matrix} have the same type, but \code{data.frame}'s
  may have mixed column types.
}
\examples{
x <- as.data.frame(cbind(1:5, LETTERS[1:5], rep(TRUE, 5)), stringsAsFactors=FALSE)
colclasses(x)
y <- as.data.frame(cbind(1:5, LETTERS[1:5], rep(TRUE, 5)), stringsAsFactors=TRUE)
colclasses(y)
}
\author{
  Mark Cowley, 18 Nov 2004
}
\seealso{
  \code{\link{'colclasses<-'}},
  \code{\link{colnameclasses}}
}
\keyword{attribute}

