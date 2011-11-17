\name{colclasses<-}
\alias{'colclasses<-'}
\alias{colclasses<-}
\title{Setter function to change the classes of the columns in x.}
\usage{
  colclasses(x) <- value
}
\arguments{
  \item{x}{a matrix or data.frame}

  \item{value}{a character vector of new column classes. if
  \code{length(value) < ncol(x)}, \code{value} is recycled}
}
\value{
  x with different classes
}
\description{
  Setter function to change the classes of the columns in
  x.
}
\examples{
y <- as.data.frame(cbind(1:5, LETTERS[1:5], rep(TRUE, 5)), stringsAsFactors=TRUE)
colclasses(y)
colclasses(y) <- c("numeric", "character", "logical")
colclasses(y)
}
\author{
  Mark Cowley
}

