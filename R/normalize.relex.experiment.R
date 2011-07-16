##' Normalize the protein ratios identified from a RelEx experiment
##' 
##' Normalize the protein ratios identified from a RelEx experiment both
##' within, and between experiments. This is a convenience function to perform the 2 most common types of
##' normalisation in the correct order.
##' More info is at ?\code{\link{normalizeWithin.relex.experiment}}
##' and ?\code{\link{normalizeBetween.relex.experiment}}.
##' 
##' @param exp A list of RelEx experiments. See \code{\link{import.relex.experiment}}.
##' @param norm.within.method The within experiment normalization method to be
##'   used. Can be one of \code{median} (default), \code{loess}, or
##'   \code{none}; see normalizeWithin.relex.experiment and
##'   normalizeWithinArrays from limma.
##' @param norm.between.method The between experiment normalization method to
##'   be used. Can be one of \code{none} (default), \code{scale}, or
##'   \code{quantile}; see \code{\link{normalizeBetween.relex.experiment}} and
##'   \code{\link[limma]{normalizeBetweenArrays}} from limma.
##' @return Returns the same list as the input, with additional normalized
##'   RelEx data and the protein level.  The list elements are as for
##'   \code{import.relex.experiment}, plus: 
##' \item{data.proteins.norm}{A data.frame of normalised protein ratios.}
##' \item{normalizations}{A list with 2 elements called \code{intra} and 
##'       \code{inter}, indicating which normalization method was used.}
##' @author Mark Cowley
##' @seealso \code{\link{normalizeWithin.relex.experiment}},
##'   \code{\link{normalizeBetween.relex.experiment}},
##'   \code{\link{import.relex.experiment}}
##' @keywords manip
##' @examples
##' # not run
##' # normalize.relex.experiment(exp, "median", "none") 
##" # normalize.relex.experiment(exp, "median", "scale")
##' # normalize.relex.experiment(exp, "loess", "scale")
##' 
normalize.relex.experiment <- function(exp, norm.within.method=c("median", "loess", "none"), norm.between.method=c("none", "scale", "quantile") ) {
	exp <- normalizeWithin.relex.experiment(exp, norm.within.method)
	exp <- normalizeBetween.relex.experiment(exp, norm.between.method)
	return( exp )
}

##' Perform intra-experiment normalization of protein ratios.
##' Perform intra-experiment normalization of protein ratios. Intra-experiment
##' normalisation is an important first step for removing skews and biases in
##' protein ratio data (see details).
##' 
##' Intra-experiment (aka within experiment) normalization is important for
##' correcting skews and biases in the protein ratio data. In a typical
##' proteomics experiment, the assumption is made that most proteins are not
##' differentially abundant, and their ratios should be centred on 1:1
##' (i.e., around the log2 ratio of 0). MA-plots are very useful for
##' highlighting skews and biases in proteomics data (see ?\code{\link{MAplot.relex}} 
##' and \code{\link{MAplot.relex.3way}}).
##' 
##' We refer to a bias as a shift away from the 1:1 ratio, and a skew as being
##' a bend in the data, such that there is a different bias dependant upon the
##' proteins' average abundance. Linear biases are best removed using median
##' normalisation, non-linear skews are best removed using loess normalization.
##' There is no need to also median normalize after loess normalization since
##' all protein ratios will be centred about the 1:1 ratio following loess.
##' 
##' This code uses \code{\link{normalizeWithinArrays}} from \code{limma},
##' 
##' @param exp A list of RelEx experiments. See \code{\link{import.relex.experiment}}
##' @param method The within experiment normalization method to be used. Can be
##'   one of \code{median} (default), \code{loess}, \code{none}; see
##'   normalizeWithinArrays from limma.
##' @param span,iterations additional arguments passed to loess. Defaults
##'   should be fine for most uses.
##' @return Returns the same list as the input, with additional normalized
##'   protein ratios.  The list elements are as for import.relex.experiment,
##'   plus: \item{data.proteins.norm}{A data.frame of normalised protein
##'   ratios.} \item{normalizations}{A list with first element called
##'   \code{intra}, indicating which normalization method was used.}
##' @author Mark Cowley
##' @seealso \code{\link{normalize.relex.experiment}},
##'   \code{\link{normalizeBetween.relex.experiment}},
##'   \code{\link{import.relex.experiment}}
##' @keywords manip
##' @examples
##' # not run
##' # normalizeWithin.relex.experiment(exp, method="median") 
##' # normalizeWithin.relex.experiment(exp, method="loess")
##' # normalizeWithin.relex.experiment(exp, method="none")
##' 
normalizeWithin.relex.experiment <- function(exp, method=c("median", "loess", "none"), span=0.3, iterations=4) {
	method <- method[1]

	if( "data.proteins.norm" %in% names(exp[[1]]) )
		warning("Normalised data will be over-written")

	require( limma )
	combo <- combine.relex.experiments(exp, normalized=FALSE)
	MA <- list(M=combo$ratio, A=combo$sn)
	MA <- new("MAList", MA)
	MAnorm <- normalizeWithinArrays(MA, method=method, span=span, iterations=iterations)
	for( i in 1:ncol(MAnorm$M) ) {
		exp[[i]]$data.proteins.norm <- exp[[i]]$data.proteins
		exp[[i]]$data.proteins.norm$ratio <- MAnorm$M[exp[[i]]$data.proteins.norm$protein, i]
		exp[[i]]$normalizations <- list(intra=method)		 
	}
	return( exp )
}

##' Perform inter-experiment normalization of protein ratios.
##' Inter-experiment normalisation should be done after intra-experiment normalization, and may
##' not be necessary. You should always be checking your data to see if a
##' particular normalization method is justified.
##' 
##' Inter-experiment (aka between experiment) normalization is important for
##' allowing comparisons between RelEx experiments which may have been
##' performed on different days, or during different runs of the MS/MS machine.
##' 
##' This should be done following intra-experiment normalisation. A boxplot of
##' the within-experiment normalised protein ratios is helpful for determining
##' the best normalisation approach. A density plot will reveal more structure
##' in the ratios; if the shapes of the ratio distributions are markedly
##' different then quantile normalisation may be in order, otherwise if the
##' shapes are all similar, but perhaps with different extents (i.e. some are
##' wider, or narrower than others) then scale normalisation may be
##' appropriate.
##' 
##' @param exp A list of RelEx experiments. See \code{\link{import.relex.experiment}}.
##' @param method The between experiment normalization method to be used. Can
##'   be one of \code{none} (default), \code{scale}, \code{quantile}; see
##'   normalizeBetweenArrays from limma.
##' @param \dots Additional arguments passed to quantile.
##' @return Returns the same list as was input, with additional normalized
##'   RelEx data and the protein level.  The list elements are as for
##'   import.relex.experiment, plus: \item{data.proteins.norm}{A data.frame of
##'   normalised protein ratios.} \item{normalizations}{A list with first
##'   elements called \code{intra}, and second element called \code{inter},
##'   indicating which normalization methods were used.}
##' @author Mark Cowley
##' @seealso \code{\link{normalize.relex.experiment}},
##'   \code{\link{normalizeWithin.relex.experiment}},
##'   \code{\link{import.relex.experiment}}
##' @keywords manip
##' @examples
##' # notrun
##' # normalizeBetween.relex.experiment(exp, method="none")
##' # normalizeBetween.relex.experiment(exp, method="scale")
##' # normalizeBetween.relex.experiment(exp, method="quantile")
##' @export
##'
normalizeBetween.relex.experiment <- function(exp, method=c("none", "scale", "quantile"), ...) {
	method <- method[1]

	if( ! "data.proteins.norm" %in% names(exp[[1]]) )
		exp <- normalizeWithin.relex.experiment(exp, method="none")

	require( limma )
	combo <- combine.relex.experiments(exp, normalized=TRUE)
	MA <- list(M=combo$ratio, A=combo$sn)
	MA <- new("MAList", MA)
	MAnorm <- normalizeBetweenArrays(MA, method=method, ...)
	for( i in 1:ncol(MAnorm$M) ) {
		exp[[i]]$data.proteins.norm$ratio <- MAnorm$M[exp[[i]]$data.proteins.norm$protein, i]
		exp[[i]]$normalizations$inter <- method
	}
	
	return( exp )
}

