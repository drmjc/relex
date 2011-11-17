\name{combine.relex.experiments}
\alias{combine.relex.experiments}
\title{Combine a list of RelEx experiments into tables for further analysis.}
\usage{
  combine.relex.experiments(x, normalized = TRUE)
}
\arguments{
  \item{x}{a list of relex experiments, produced by
  \code{\link{import.relex.experiment}}. Each list element
  corresponds to an individual RelEx experiment, and is a
  list containing at least the table \code{data.proteins}
  or \code{data.proteins.norm}.}

  \item{normalized}{logical: If TRUE, then the
  \code{data.proteins.norm} table from each experiment is
  used, else the \code{data.proteins} table is used.
  Defaults to TRUE.}
}
\value{
  a list with the following elements: \item{ratio}{a
  data.frame of protein ratios (log2 scale). Usually this
  will be the normalised ratios.} \item{sn}{a data.frame of
  Signal/Noise ratios (log10 scale)} \item{ratioSD}{a
  data.frame of the standard deviation of the normalised
  protein ratios (obtained from the peptides)}
  \item{snSD}{a data.frame of the standard deviation of the
  Signal/Noise ratios (obtained from the peptides)}
  \item{proteins}{A vector of protein names that were
  identified in any experiment} \item{protein.counts}{A
  named numeric vector of counts, indicating in how many
  experiments each protein was observered, where the names
  are the protein names.}
}
\description{
  Combine a list of relex experiments into 4 tables:
  \dQuote{ratio}, \dQuote{sn} (Signal to Noise),
  \dQuote{ratioSD} (stddev of the ratios), \dQuote{snSD}
  (stdev of the Signal/Noise ratios), and 2 vectors:
  \dQuote{proteins} and \dQuote{protein.counts}. The
  standard deviations are calculated from the peptides that
  form each protein.
}
\examples{
\dontrun{
combine.relex.experiments(x)
combine.relex.experiments(x, normalized = FALSE)
}
}
\author{
  Mark Cowley
}
\seealso{
  \code{\link{import.relex.experiment}}
}
\keyword{manip}

