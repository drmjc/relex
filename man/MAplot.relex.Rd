\name{MAplot.relex}
\alias{MAplot.relex}
\title{Perform a ratio-vs-average plot of peptide or protein ratios.
Plot a ratio-vs-average (MA-plot) of peptide ratios, or protein ratios, or
normalised protein ratios vs the log10 signal/noise ratio (an estimate of
the amount of expression). \code{\link{MAplot.relex.3way}} plots all three.}
\usage{
  MAplot.relex(x, what = c("proteins.norm", "peptides",
  "proteins"), main = "", ylim = NULL, loess = TRUE, legend
  = TRUE, ...)
}
\arguments{
  \item{x}{A list of RelEx experiments.}

  \item{what}{Which data to plot: \code{proteins.norm} (the
  default), \code{peptides}, or \code{proteins}.}

  \item{main}{Arguments passed on to to
  \code{\link{plot}}.}

  \item{ylim}{If NULL then the ylim is determined
  automatically, otherwise it can be over-ridden here.}

  \item{loess}{logical: add a loess smothed line over the
  MA plot in purple?}

  \item{legend}{logical: add a legend to the plot,
  indicating the number of proteins or peptides in the
  plot?}

  \item{\dots}{other arguments passed to plot}
}
\value{
  none. generates an MA plot
}
\description{
  MA-plots are extremely useful for highlighting skews and
  biases in a comparison of 2 biological samples, such as 2
  samples labelled with N14 vs N15. The user may select
  which quantity to plot (peptides/proteins/normalised
  proteins), and a loess curve in pink indicates any
  systematic trends that may exist in the data.
}
\examples{
\dontrun{
MAplot.relex(exp, "peptides")
MAplot.relex(exp, "proteins") MAplot.relex(exp, "proteins.norm")
MAplot.relex(exp, "proteins.norm", loess=FALSE, legend=FALSE)
MAplot.relex.3way(exp) MAplot.relex.3way(exp, ylim = c(-4, 4))
}
}
\author{
  Mark Cowley
}
\seealso{
  \code{\link{MAplot.relex.3way}}, \code{\link{lowess}},
  \code{\link{symmetricise}}
}
\keyword{hplot}

