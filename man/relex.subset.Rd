\name{relex.subset}
\alias{relex.subset}
\title{Choose a subset of samples from a set of relex experiments...}
\usage{relex.subset(x, cols, rmNArows=TRUE)
}
\description{Choose a subset of samples from a set of relex experiments
Choose a subset of samples from a set of relex experiments. Useful if you
decide a sample is an outlier, and want to remove it, or if you want to do
some sort of leave-1-out analysis.}
\value{a list much like the input though with fewer columns.}
\author{Mark Cowley}
\seealso{\code{\link{rmNArows.relex}}}
\keyword{manip}
\arguments{\item{x}{A list of RelEx experiments.}
\item{cols}{A vector of column indices, or names of the columns to keep.}
\item{rmNArows}{when you subset columns, you can sometimes end up with rows
that are all NA. setting rmNArows to TRUE removes those rows, FALSE means
leave them in. See rmNArows.relex}
}
\examples{# not run
# relex.subset(x, c(1,3,4,5), rmNArows = TRUE)
# relex.subset(x, c("A1", "B1", "C1", "D1", "E1"), rmNArows = TRUE)}

