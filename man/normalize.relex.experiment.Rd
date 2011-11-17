\name{normalize.relex.experiment}
\alias{normalize.relex.experiment}
\title{Normalize the protein ratios identified from a RelEx experiment}
\usage{
  normalize.relex.experiment(exp, norm.within.method =
  c("median", "loess", "none"), norm.between.method =
  c("none", "scale", "quantile"))
}
\arguments{
  \item{exp}{A list of RelEx experiments. See
  \code{\link{import.relex.experiment}}.}

  \item{norm.within.method}{The within experiment
  normalization method to be used. Can be one of
  \code{median} (default), \code{loess}, or \code{none};
  see normalizeWithin.relex.experiment and
  normalizeWithinArrays from limma.}

  \item{norm.between.method}{The between experiment
  normalization method to be used. Can be one of
  \code{none} (default), \code{scale}, or \code{quantile};
  see \code{\link{normalizeBetween.relex.experiment}} and
  \code{\link[limma]{normalizeBetweenArrays}} from limma.}
}
\value{
  Returns the same list as the input, with additional
  normalized RelEx data and the protein level.  The list
  elements are as for \code{import.relex.experiment}, plus:
  \item{data.proteins.norm}{A data.frame of normalised
  protein ratios.} \item{normalizations}{A list with 2
  elements called \code{intra} and \code{inter}, indicating
  which normalization method was used.}
}
\description{
  Normalize the protein ratios identified from a RelEx
  experiment both within, and between experiments. This is
  a convenience function to perform the 2 most common types
  of normalisation in the correct order. More info is at
  ?\code{\link{normalizeWithin.relex.experiment}} and
  ?\code{\link{normalizeBetween.relex.experiment}}.
}
\examples{
\dontrun{
normalize.relex.experiment(exp, "median", "none")
normalize.relex.experiment(exp, "median", "scale")
normalize.relex.experiment(exp, "loess", "scale")
}
}
\author{
  Mark Cowley
}
\seealso{
  \code{\link{normalizeWithin.relex.experiment}},
  \code{\link{normalizeBetween.relex.experiment}},
  \code{\link{import.relex.experiment}}
}
\keyword{manip}

