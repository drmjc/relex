\name{relex.merge}
\alias{relex.merge}
\title{Merge multiple relex experiments together
Once multiple csv files from the same group of experiments have been
combined using \code{\link{combine.relex.experiments}}, you can merge these
multiple experiments.}
\usage{
  relex.merge(x, y)
}
\arguments{
  \item{x,y}{each is the output from
  \code{combine.relex.experiments}, which is a list
  containing at least the tables sn, ratio, snSD, ratioSD.
  Each table's rownames must be the protein names, and the
  colnames must be the experiment name.}
}
\value{
  A list with the same elements as the input lists, but the
  tables have been cbinded together (x then y), and the
  proteins and protein.counts are updated.
}
\description{
  Merge multiple relex experiments together Once multiple
  csv files from the same group of experiments have been
  combined using \code{\link{combine.relex.experiments}},
  you can merge these multiple experiments.
}
\examples{
\dontrun{
data.merged <- relex.merge(data.10vs30, data.30vs10)
}
}
\author{
  Mark Cowley
}
\seealso{
  \code{\link{combine.relex.experiments}}
}
\keyword{manip}

