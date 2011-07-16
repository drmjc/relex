##' Import one or more csv file(s) produced by RelEx
##' 
##' This is the primary function to import RelEx data.  You can import multiple files at
##' once. Check the example section for examples of these RelEx formatted csv files.
##' 
##' @param filename a character vector of filenames. These should include the
##'   full path.
##' @param normalize if TRUE, then a median intra-array, and no inter-array
##'   normalization is performed. If you'd like any other normalization
##'   options, set this to FALSE, then subseqently normalize using
##'   normalize.relex.experiment. Defaults to FALSE.
##' @param verbose if TRUE then informative messages are printed. Defaults to
##'   TRUE.
##' @return A list of relex experiments are returned, 1 element per file.  Each
##'   experiment is a list of 3-4 objects: \item{data.peptides}{a data.frame of
##'   ratios, one row per peptide.} \item{data.proteins}{a data.frame of
##'   ratios, one row per protein} \item{proteins}{A character vector of
##'   protein names.} \item{data.protein.norm}{[optional] a data.frame of
##'   normalised ratios, one row per protein.} \item{normalizations}{[optional]
##'   a list indicating which normalizations have been applied.}
##' @note If this function does not work with your type of input file, please
##'   contact the package maintainer.  Better yet, if you can write an importer
##'   for your particular type of file, then also contact the package
##'   maintainer.
##' @author Mark Cowley, 11/7/07
##' @seealso \code{\link{normalize.relex.experiment}},
##'   \code{\link{normalizeWithin.relex.experiment}},
##'   \code{\link{normalizeBetween.relex.experiment}},
##'   \code{\link{combine.relex.experiments}}
##' @keywords IO file
##' @examples
##' # not run
##' # import.relex.experiment(filename)
##' # import.relex.experiment(filename, normalize = TRUE)
##' # import.relex.experiment(filename, normalize = FALSE, verbose = FALSE)
##' # import.relex.experiment(dir("/path/to/files", pattern=".*csv"), normalize = FALSE, verbose = TRUE)
##' @export
##'
import.relex.experiment <- function(filename, normalize=FALSE, verbose=TRUE) {
	stopifnot( all(file.exists(filename)) )

	if( length(filename) > 1 ) {
		res <- list()
		for(i in 1:length(filename)) {
			if( verbose ) cat("importing", filename[i], "\n")
			res[[i]] <- import.relex.experiment(filename[i], normalize=normalize, verbose=verbose)
		}
		names(res) <- strip.fileextension( basename(filename) )
		return( res )
	}
	else {
		if( verbose ) cat("importing\n")
		map <- read.delim(filename, as.is=TRUE, header=FALSE)

		data.peptides <- as.data.frame(matrix(0, nrow=sum(map[,1] %in% c('S', 'R')), ncol=4), stringsAsFactors=FALSE) # the rows start with S or R for Sample, Reference.
		colnames(data.peptides) <- c("protein", "peptide", "ratio", "noise")
		colclasses(data.peptides) <- c("character", "character", "numeric", "numeric")

		if( verbose ) cat("retrieving peptide ratios\n")
		protein <- ""
		peptide <- ""
		idx <- 1
		for(i in 1:nrow(map)) {
			if(map[i,1] == "P")
				protein <- map[i,2]
			else if( map[i,1] %in% c("S", "R") ) {
				data.peptides[idx,] <- c(protein, map[i,3], map[i,4], map[i,7])
				idx <- idx + 1
			}
		}
		colclasses(data.peptides) <- c("character", "character", "numeric", "numeric")

		if( verbose ) cat("summarising peptides --> proteins\n")
		#
		# convert peptides into protein measurements
		#
		data.peptides$ratio <- log2(data.peptides$ratio)
		data.peptides$noise <- log10(data.peptides$noise)

		data.proteins <- aggregate(data.peptides[,c(3,4)], list(data.peptides$protein), mean)
		colnames(data.proteins)[1] <- "protein"
		data.proteins[,1] <- as.character(data.proteins[,1])

		tmp <- aggregate(data.peptides[,c(3,4)], list(data.peptides$protein), sd)[,2:3]
		colnames(tmp) <- c("ratioSD", "noiseSD")
		data.proteins <- cbind(data.proteins, tmp)


		exp <- list(data.peptides=data.peptides,
					 data.proteins=data.proteins,
					 proteins=data.proteins[,1])

		if( normalize ) {
			if( verbose ) cat("normalising protein ratios, using default settings\n")
			exp <- normalize.relex.experiment(exp)
		}

		return( exp )
	}
}
