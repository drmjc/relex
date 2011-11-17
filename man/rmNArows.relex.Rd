\name{rmNArows.relex}
\alias{rmNArows.relex}
\title{Removes rows from a RelEx experiment that have zero detection counts.
Removes rows from a RelEx experiment that have zero detection counts. This
will update all tables in the list, as well as the proteins and
protein.counts vectors}
\usage{
  rmNArows.relex(x)
}
\arguments{
  \item{x}{A list of RelEx experiments.}
}
\value{
  a list much like \code{x} but where
  \code{min(protein.counts) >= 1}
}
\description{
  Removes rows from a RelEx experiment that have zero
  detection counts. Removes rows from a RelEx experiment
  that have zero detection counts. This will update all
  tables in the list, as well as the proteins and
  protein.counts vectors
}
\author{
  Mark Cowley
}
\keyword{manip}

