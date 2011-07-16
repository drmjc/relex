##' Export a Relex experiment to tab delimited files.
##'
##' Export the peptide, protein and normalised protein ratios from a Relex
##' experiment as tab delimited files.
##' 
##' @param exp A list of RelEx experiments
##' @param dir The path to the directory where you'd like the files to be made.
##' @param prefix The file name's prefix. You can leave this as \dQuote{}, or
##'   name it something useful that will separate these files that will be made
##'   from others that you may wish to make.
##' @return this exports 2 or 3 tab delimited files to a directory,
##'   corresponding to the peptides, proteins and normalised protein ratios.
##' @author Mark Cowley
##' @seealso \code{\link{write.delim}}
##' @keywords IO file
##' @examples
##' # not run
##' # export.relex.experiment(exp, dir, prefix) 
##' # export.relex.experiment(exp, "/path/to/output", "")
##' # export.relex.experiment(exp, "/path/to/output", "MUTvsWT")
##' @export
##'
export.relex.experiment <- function(exp, dir="./", prefix="") {
	f <- file.path(dir, paste(prefix, c("peptides.csv", "proteins.csv", "proteins.norm.csv"), sep="."))
	write.delim(exp$data.peptides, f[1])
	write.delim(exp$data.proteins, f[2])
	if( "data.proteins.norm" %in% names(exp) )
		write.delim(exp$data.proteins.norm, f[3])
}
