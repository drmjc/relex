##' Perform a ratio-vs-average plot of peptide or protein ratios.
##' Plot a ratio-vs-average (MA-plot) of peptide ratios, or protein ratios, or
##' normalised protein ratios vs the log10 signal/noise ratio (an estimate of
##' the amount of expression). \code{\link{MAplot.relex.3way}} plots all three.
##' 
##' MA-plots are extremely useful for highlighting skews and biases in a
##' comparison of 2 biological samples, such as 2 samples labelled with N14 vs
##' N15. The user may select which quantity to plot
##' (peptides/proteins/normalised proteins), and a loess curve in pink
##' indicates any systematic trends that may exist in the data.
##' 
##' @param x A list of RelEx experiments.
##' @param what Which data to plot: \code{proteins.norm} (the default),
##'   \code{peptides}, or \code{proteins}.
##' @param main Arguments passed on to to \code{\link{plot}}.
##' @param ylim If NULL then the ylim is determined automatically, otherwise it
##'   can be over-ridden here.
##' @param loess logical: add a loess smothed line over the MA plot in purple?
##' @param legend logical: add a legend to the plot, indicating the number of proteins or peptides in the plot?
##' @param \dots other arguments passed to plot
##' @return none. generates an MA plot
##' @author Mark Cowley
##' @seealso \code{\link{MAplot.relex.3way}}, \code{\link{lowess}}, \code{\link{symmetricise}}
##' @keywords hplot
##' @examples
##' # not run
##' # MAplot.relex(exp, "peptides")
##' # MAplot.relex(exp, "proteins") MAplot.relex(exp, "proteins.norm")
##' # MAplot.relex(exp, "proteins.norm", loess=FALSE, legend=FALSE)
##' # MAplot.relex.3way(exp) MAplot.relex.3way(exp, ylim = c(-4, 4))
##' @export
MAplot.relex <- function(x, what=c("proteins.norm", "peptides", "proteins"), main="", ylim=NULL, loess=TRUE, legend=TRUE, ...) {
	if( what == "peptides") x <- x$data.peptides
	else if( what == "proteins") x <- x$data.proteins
	else if( what == "proteins.norm") x <- x$data.proteins.norm
	else stop("You must choose one of 'peptides', 'proteins', or 'proteins.norm'.\n")
	
	if( is.null(ylim) )
		ylim <- symmetricise(range(x$ratio))

	plot(x$noise, x$ratio, xlab="S/N (log10)", ylab="Relex ratio (log2)", main=main, ylim=ylim, ...)
	if( loess ) lines(lowess(x$noise, x$ratio), col=6)
	if( legend ) legend("bottomright", paste("n =", nrow(x)), bty="n")
}


##' MA-plot of relex data in three panels.
##' Plot a ratio-vs-average (MA-plot) of peptide ratios, protein ratios, and
##' normalised protein ratios vs the log10 signal/noise ratio (an estimate of
##' the amount of expression).
##' 
##' MA-plots are extremely useful for highlighting skews and biases in a
##' comparison of 2 biological samples, such as 2 samples labelled with N14 vs
##' N15. The user may select which quantity to plot
##' (peptides/proteins/normalised proteins), and a loess curve in pink
##' indicates any systematic trends that may exist in the data.
##' 
##' @param x A list of RelEx experiments.
##' @param ylim If NULL then the ylim is determined automatically, otherwise it
##'   can be over-ridden here.
##' @param loess logical: add a loess smothed line over the MA plot in purple?
##' @param legend logical: add a legend to the plot, indicating the number of proteins or peptides in the plot?
##' @param do.par logical: set the \code{par(mfrow)} to \code{c(1,3)}?
##' @param \dots other arguments passed to plot
##' @return none. generates three MA plots
##' @author Mark Cowley
##' @seealso \code{\link{MAplot.relex}}, \code{\link{lowess}}, \code{\link{symmetricise}}
##' @keywords hplot
##' @examples
##' # not run
##' # MAplot.relex(exp, "peptides")
##' # MAplot.relex(exp, "proteins") MAplot.relex(exp, "proteins.norm")
##' # MAplot.relex(exp, "proteins.norm", loess=FALSE, legend=FALSE)
##' # MAplot.relex.3way(exp) MAplot.relex.3way(exp, ylim = c(-4, 4))
##' @export
MAplot.relex.3way <- function(x, ylim=c(-4,4), loess=TRUE, legend=TRUE, do.par=TRUE, ...) {
	if( do.par ) par(mfrow=c(1,3))
	MAplot.relex(x, what="peptides", main="Peptide ratios", ylim=ylim, loess=loess, legend=legend, ...)
	MAplot.relex(x, what="proteins", main="Protein ratios", ylim=ylim, loess=loess, legend=legend, ...)
	MAplot.relex(x, what="proteins.norm", main="Normalised Protein ratios", ylim=ylim, loess=loess, legend=legend, ...)
}
