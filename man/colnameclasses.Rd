\name{colnameclasses}
\alias{colnameclasses}
\title{Print out the colnames and their classes in a 2 column table...}
\usage{colnameclasses(x)
}
\description{Print out the colnames and their classes in a 2 column table}
\value{a 2 column table of colname, and colclasses}
\author{Mark Cowley, 27/3/07}
\seealso{\code{\link{colclasses}}}
\arguments{\item{x}{a table}
}
\examples{y <- as.data.frame(cbind(1:5, LETTERS[1:5], rep(TRUE, 5)), stringsAsFactors=TRUE)
colnameclasses(y)}

