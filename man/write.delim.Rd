\name{write.delim}
\alias{write.delim}
\title{Write a matrix-like object as a tab-delimited txt file.}
\usage{write.delim(x, file, na="---", row.names=FALSE, col.names=TRUE, append=FALSE, ...)
}
\description{Write a matrix-like object as a tab-delimited txt file.
\code{write.delim} writes tab delimited files, but offers some enhancements
relative to the standard read.table. These improvements mostly relate to
enhanced control over the row names that are exported.}
\details{If you use the default options, then it works nicely with
\code{\link{read.delim}}.}
\value{a tab-delimited txt file is written out.}
\author{Mark Cowley}
\seealso{\code{\link{write.table}}, \code{\link{read.delim}}}
\keyword{IO}
\keyword{file}
\arguments{\item{x}{a \code{matrix} or \code{data.frame}}
\item{file}{The file name (including path) which will be created. It will be
over-written without warning.}
\item{na}{What to replace NA's with. Sensible options include \dQuote{---}
(the default), \dQuote{}, or \dQuote{-}.}
\item{row.names}{either TRUE, FALSE, NULL, NA, or a character(1).  if FALSE,
then no row names will be written out, if TRUE then row.names will be
written out, with a blank column header (i.e. dQuote); if NULL or NA, then
row.names will be written out with no column name (just like
write.table(row.names=TRUE) would do) if row.names is a character vector
of length 1, then write the row.names out, and use use the supplied word
as the column name for the rownames else error}
\item{col.names}{TRUE or FALSE. If TRUE then write out the column names.}
\item{append}{append to end, or overwrite the file?}
\item{\dots}{Further arguments passed to \code{\link{write.table}}}
}

