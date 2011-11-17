\name{normalizeBetween.relex.experiment}
\alias{normalizeBetween.relex.experiment}
\title{Perform inter-experiment normalization of protein ratios.
Inter-experiment normalisation should be done after intra-experiment normalization, and may
not be necessary. You should always be checking your data to see if a
particular normalization method is justified.}
\usage{
  normalizeBetween.relex.experiment(exp, method = c("none",
  "scale", "quantile"), ...)
}
\arguments{
  \item{exp}{A list of RelEx experiments. See
  \code{\link{import.relex.experiment}}.}

  \item{method}{The between experiment normalization method
  to be used. Can be one of \code{none} (default),
  \code{scale}, \code{quantile}; see normalizeBetweenArrays
  from limma.}

  \item{\dots}{Additional arguments passed to quantile.}
}
\value{
  Returns the same list as was input, with additional
  normalized RelEx data and the protein level.  The list
  elements are as for import.relex.experiment, plus:
  \item{data.proteins.norm}{A data.frame of normalised
  protein ratios.} \item{normalizations}{A list with first
  elements called \code{intra}, and second element called
  \code{inter}, indicating which normalization methods were
  used.}
}
\description{
  Inter-experiment (aka between experiment) normalization
  is important for allowing comparisons between RelEx
  experiments which may have been performed on different
  days, or during different runs of the MS/MS machine.
}
\details{
  This should be done following intra-experiment
  normalisation. A boxplot of the within-experiment
  normalised protein ratios is helpful for determining the
  best normalisation approach. A density plot will reveal
  more structure in the ratios; if the shapes of the ratio
  distributions are markedly different then quantile
  normalisation may be in order, otherwise if the shapes
  are all similar, but perhaps with different extents (i.e.
  some are wider, or narrower than others) then scale
  normalisation may be appropriate.
}
\examples{
\dontrun{
normalizeBetween.relex.experiment(exp, method="none")
normalizeBetween.relex.experiment(exp, method="scale")
normalizeBetween.relex.experiment(exp, method="quantile")
}
}
\author{
  Mark Cowley
}
\seealso{
  \code{\link{normalize.relex.experiment}},
  \code{\link{normalizeWithin.relex.experiment}},
  \code{\link{import.relex.experiment}}
}
\keyword{manip}

