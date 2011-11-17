#' Removes rows from a RelEx experiment that have zero detection counts.
#' Removes rows from a RelEx experiment that have zero detection counts. This
#' will update all tables in the list, as well as the proteins and
#' protein.counts vectors
#' 
#' @param x A list of RelEx experiments.
#' @return a list much like \code{x} but where \code{min(protein.counts) >= 1}
#' @author Mark Cowley
#' @keywords manip
#' @export
#'
rmNArows.relex <- function(x) {
	idx <- which(x$protein.counts > 0)
	if( length(idx) < length(x$protein.counts) ) {
		x$sn <- x$sn[idx,]
		x$ratio <- x$ratio[idx,]
		x$snSD <- x$snSD[idx,]
		x$ratioSD <- x$ratioSD[idx,]
		x$proteins <- x$proteins[idx]
		x$protein.counts <- x$protein.counts[idx]
	}
	return( x )
}
