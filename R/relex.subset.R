#' Choose a subset of samples from a set of relex experiments
#' Choose a subset of samples from a set of relex experiments. Useful if you
#' decide a sample is an outlier, and want to remove it, or if you want to do
#' some sort of leave-1-out analysis.
#' 
#' @param x A list of RelEx experiments.
#' @param cols A vector of column indices, or names of the columns to keep.
#' @param rmNArows when you subset columns, you can sometimes end up with rows
#'   that are all NA. setting rmNArows to TRUE removes those rows, FALSE means
#'   leave them in. See rmNArows.relex
#' @return a list much like the input though with fewer columns.
#' @author Mark Cowley
#' @seealso \code{\link{rmNArows.relex}}
#' @keywords manip
#' @examples
#' \dontrun{
#' relex.subset(x, c(1,3,4,5), rmNArows = TRUE)
#' relex.subset(x, c("A1", "B1", "C1", "D1", "E1"), rmNArows = TRUE)
#' }
#' @export
#'
relex.subset <- function(x, cols, rmNArows=TRUE) {
	if( missing(cols) )
		stop("You must choose some columns to retain.\n")
	res <- x
	res$sn <- res$sn[,cols]
	res$ratio <- res$ratio[,cols]
	res$snSD <- res$snSD[,cols]
	res$ratioSD <- res$ratioSD[,cols]
	res$protein.counts <- apply(!is.na(res$ratio), 1, sum)
	# res$proteins remains unchanged

	if( rmNArows && any(res$protein.counts == 0) ) {
		res <- rmNArows.relex(res)
	}
	return( res )
}
