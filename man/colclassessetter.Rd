\name{colclassessetter}
\alias{colclasses<-}
\title{Setter function to change the classes of the columns in x.}
\usage{colclasses(x) <- value}
\description{Setter function to change the classes of the columns in x.}
\value{x with different classes}
\author{Mark Cowley}
\alias{'colclasses<-'}
\arguments{\item{x}{a matrix or data.frame}
\item{value}{a character vector of new column classes. if \code{length(value) < ncol(x)}, \code{value} is recycled}
}
\examples{y <- as.data.frame(cbind(1:5, LETTERS[1:5], rep(TRUE, 5)), stringsAsFactors=TRUE)
colclasses(y)
colclasses(y) <- c("numeric", "character", "logical")
colclasses(y)}

