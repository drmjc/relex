##' Combine a list of RelEx experiments into tables for further analysis.
##'
##' Combine a list of relex experiments into 4 tables: \dQuote{ratio}, \dQuote{sn} (Signal to
##' Noise), \dQuote{ratioSD} (stddev of the ratios), \dQuote{snSD} (stdev of the Signal/Noise
##' ratios), and 2 vectors: \dQuote{proteins} and \dQuote{protein.counts}. The standard
##' deviations are calculated from the peptides that form each protein.
##' 
##' @param x a list of relex experiments, produced by \code{\link{import.relex.experiment}}.
##'   Each list element corresponds to an individual RelEx experiment, and is a
##'   list containing at least the table \code{data.proteins} or
##'   \code{data.proteins.norm}.
##' @param normalized logical: If TRUE, then the \code{data.proteins.norm} table from
##'   each experiment is used, else the \code{data.proteins} table is used.
##'   Defaults to TRUE.
##' @return a list with the following elements: \item{ratio}{a data.frame of
##'   protein ratios (log2 scale). Usually this will be the normalised ratios.}
##'   \item{sn}{a data.frame of Signal/Noise ratios (log10 scale)}
##'   \item{ratioSD}{a data.frame of the standard deviation of the normalised
##'   protein ratios (obtained from the peptides)} \item{snSD}{a data.frame of
##'   the standard deviation of the Signal/Noise ratios (obtained from the
##'   peptides)} \item{proteins}{A vector of protein names that were identified
##'   in any experiment} \item{protein.counts}{A named numeric vector of
##'   counts, indicating in how many experiments each protein was observered,
##'   where the names are the protein names.}
##' @author Mark Cowley
##' @seealso \code{\link{import.relex.experiment}}
##' @keywords manip
##' @examples
##' #not run
##' # combine.relex.experiments(x)
##' # combine.relex.experiments(x, normalized = FALSE)
##' @export
##'
combine.relex.experiments <- function(x, normalized=TRUE) {
	if( normalized && "data.proteins.norm" %in% names(x[[1]]))
		idx <- match("data.proteins.norm", names(x[[1]]))
	else
		idx <- match("data.proteins", names(x[[1]]))

	prot.per.exp <- lapply(x, function(y) y[[idx]]$protein)
	unique.proteins <- sort( unionN( prot.per.exp ) )

	tmp <- matrix(NA, nrow=length(unique.proteins), ncol=length(x))
	rownames(tmp) <- unique.proteins
	colnames(tmp) <- names(x)
	tmp <- as.data.frame(tmp, stringsAsFactors=FALSE)

	sn <- tmp
	ratio <- tmp
	snSD <- tmp
	ratioSD <- tmp

	for(i in 1:length(x)) {
		sn[x[[i]][[idx]]$protein, i] <- x[[i]][[idx]]$noise
		ratio[x[[i]][[idx]]$protein, i] <- x[[i]][[idx]]$ratio
		snSD[x[[i]][[idx]]$protein, i] <- x[[i]][[idx]]$noiseSD
		ratioSD[x[[i]][[idx]]$protein, i] <- x[[i]][[idx]]$ratioSD
	}

	p.counts <- ucounts( unlist(prot.per.exp) )
	p.counts <- p.counts[order(names(p.counts))]

	res <- list(sn=sn, ratio=ratio, snSD=snSD, ratioSD=ratioSD, proteins=unique.proteins, protein.counts=p.counts)

	return( res )
}
