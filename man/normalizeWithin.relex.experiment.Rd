\name{normalizeWithin.relex.experiment}
\alias{normalizeWithin.relex.experiment}
\title{Perform intra-experiment normalization of protein ratios.}
\usage{normalizeWithin.relex.experiment(exp, method=c("median", "loess", "none"), span=0.3, iterations=4)
}
\description{Perform intra-experiment normalization of protein ratios.
Perform intra-experiment normalization of protein ratios. Intra-experiment
normalisation is an important first step for removing skews and biases in
protein ratio data (see details).}
\details{Intra-experiment (aka within experiment) normalization is important for
correcting skews and biases in the protein ratio data. In a typical
proteomics experiment, the assumption is made that most proteins are not
differentially abundant, and their ratios should be centred on 1:1
(i.e., around the log2 ratio of 0). MA-plots are very useful for
highlighting skews and biases in proteomics data (see ?\code{\link{MAplot.relex}} 
and \code{\link{MAplot.relex.3way}}).

We refer to a bias as a shift away from the 1:1 ratio, and a skew as being
a bend in the data, such that there is a different bias dependant upon the
proteins' average abundance. Linear biases are best removed using median
normalisation, non-linear skews are best removed using loess normalization.
There is no need to also median normalize after loess normalization since
all protein ratios will be centred about the 1:1 ratio following loess.

This code uses \code{\link{normalizeWithinArrays}} from \code{limma},}
\value{Returns the same list as the input, with additional normalized
protein ratios.  The list elements are as for import.relex.experiment,
plus: \item{data.proteins.norm}{A data.frame of normalised protein
ratios.} \item{normalizations}{A list with first element called
\code{intra}, indicating which normalization method was used.}}
\author{Mark Cowley}
\seealso{\code{\link{normalize.relex.experiment}},
\code{\link{normalizeBetween.relex.experiment}},
\code{\link{import.relex.experiment}}}
\keyword{manip}
\arguments{\item{exp}{A list of RelEx experiments. See \code{\link{import.relex.experiment}}}
\item{method}{The within experiment normalization method to be used. Can be
one of \code{median} (default), \code{loess}, \code{none}; see
normalizeWithinArrays from limma.}
\item{span,iterations}{additional arguments passed to loess. Defaults
should be fine for most uses.}
}
\examples{# not run
# normalizeWithin.relex.experiment(exp, method="median") 
# normalizeWithin.relex.experiment(exp, method="loess")
# normalizeWithin.relex.experiment(exp, method="none")}

