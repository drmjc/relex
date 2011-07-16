##' Merge multiple relex experiments together
##' Once multiple csv files from the same group of experiments have been
##' combined using \code{\link{combine.relex.experiments}}, you can merge these
##' multiple experiments.
##' 
##' @param x,y each is the output from \code{combine.relex.experiments}, which
##'   is a list containing at least the tables sn, ratio, snSD, ratioSD. Each
##'   table's rownames must be the protein names, and the colnames must be the
##'   experiment name.
##' @return A list with the same elements as the input lists, but the tables
##'   have been cbinded together (x then y), and the proteins and
##'   protein.counts are updated.
##' @author Mark Cowley
##' @seealso \code{\link{combine.relex.experiments}}
##' @keywords manip
##' @examples
##' # not run
##' # data.merged <- relex.merge(data.10vs30, data.30vs10)
##' @export
relex.merge <- function(x, y) {
	uprots <- sort( union(rownames(x$proteins), rownames(y$proteins)) )

	res <- list()
	res$sn <- cbind.smart(x$sn, y$sn)
	res$ratio <- cbind.smart(x$ratio, y$ratio)
	res$snSD <- cbind.smart(x$snSD, y$snSD)
	res$ratioSD <- cbind.smart(x$ratioSD, y$ratioSD)

	res$proteins <- rownames(res$sn)
	res$protein.counts <- apply(res$ratio, 1, function(x) ncol(res$ratio) - sum(is.na(x)))

	return( res )
}
