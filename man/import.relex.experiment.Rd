\name{import.relex.experiment}
\alias{import.relex.experiment}
\title{Import one or more csv file(s) produced by RelEx}
\usage{
  import.relex.experiment(filename, normalize = FALSE,
  verbose = TRUE)
}
\arguments{
  \item{filename}{a character vector of filenames. These
  should include the full path.}

  \item{normalize}{if TRUE, then a median intra-array, and
  no inter-array normalization is performed. If you'd like
  any other normalization options, set this to FALSE, then
  subseqently normalize using normalize.relex.experiment.
  Defaults to FALSE.}

  \item{verbose}{if TRUE then informative messages are
  printed. Defaults to TRUE.}
}
\value{
  A list of relex experiments are returned, 1 element per
  file.  Each experiment is a list of 3-4 objects:
  \item{data.peptides}{a data.frame of ratios, one row per
  peptide.} \item{data.proteins}{a data.frame of ratios,
  one row per protein} \item{proteins}{A character vector
  of protein names.} \item{data.protein.norm}{[optional] a
  data.frame of normalised ratios, one row per protein.}
  \item{normalizations}{[optional] a list indicating which
  normalizations have been applied.}
}
\description{
  This is the primary function to import RelEx data.  You
  can import multiple files at once. Check the example
  section for examples of these RelEx formatted csv files.
}
\note{
  If this function does not work with your type of input
  file, please contact the package maintainer.  Better yet,
  if you can write an importer for your particular type of
  file, then also contact the package maintainer.
}
\examples{
pkg <- file.path(path.package("relex"), "data")
filename <- dir(pkg, pattern="30-30", full=TRUE)[1]
import.relex.experiment(filename)
\dontrun{
import.relex.experiment(filename, normalize = TRUE)
import.relex.experiment(filename, normalize = FALSE, verbose = FALSE)
import.relex.experiment(dir("/path/to/files", pattern=".*csv"), normalize = FALSE, verbose = TRUE)
}
}
\author{
  Mark Cowley, 11/7/07
}
\seealso{
  \code{\link{normalize.relex.experiment}},
  \code{\link{normalizeWithin.relex.experiment}},
  \code{\link{normalizeBetween.relex.experiment}},
  \code{\link{combine.relex.experiments}}
}
\keyword{file}
\keyword{IO}

