\name{export.relex.experiment}
\alias{export.relex.experiment}
\title{Export a Relex experiment to tab delimited files.}
\usage{
  export.relex.experiment(exp, dir = "./", prefix = "")
}
\arguments{
  \item{exp}{A list of RelEx experiments}

  \item{dir}{The path to the directory where you'd like the
  files to be made.}

  \item{prefix}{The file name's prefix. You can leave this
  as \dQuote{}, or name it something useful that will
  separate these files that will be made from others that
  you may wish to make.}
}
\value{
  this exports 2 or 3 tab delimited files to a directory,
  corresponding to the peptides, proteins and normalised
  protein ratios.
}
\description{
  Export the peptide, protein and normalised protein ratios
  from a Relex experiment as tab delimited files.
}
\examples{
\dontrun{
export.relex.experiment(exp, dir, prefix)
export.relex.experiment(exp, "/path/to/output", "")
export.relex.experiment(exp, "/path/to/output", "MUTvsWT")
}
}
\author{
  Mark Cowley
}
\seealso{
  \code{\link{write.delim}}
}
\keyword{file}
\keyword{IO}

