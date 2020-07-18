\name{readExpData}
\Rdversion{1.1}
\alias{readExpData}

\title{
Read in expression data
}
\description{
This is a wrapper function of \code{read.delim} for reading in
expression data matrix in tab-delimited format.
}
\usage{
readExpData(file = "arrayData.txt", ...)
}

\arguments{
  \item{file}{
character string, the full path name to the expression data file in
tab-delimited format. Rows are genes, columns are array samples.
}
\item{\dots}{
    other arguments to be passed into \code{read.delim} function.
  }
}
\details{
\code{readExpData} is a wrapper function of \code{read.delim}. Please
check help information of \code{read.delim} for more details.
}
\value{
  A data.frame (matrix-like) of gene expression data. Rows are genes,
  columns are array samples.
}
\references{
  Luo, W., Friedman, M., Shedden K., Hankenson, K. and Woolf, P GAGE:
  Generally Applicable Gene Set Enrichment for Pathways Analysis. BMC
  Bioinformatics 2009, 10:161
}
\author{
  Weijun Luo <luo_weijun@yahoo.com>
}

\seealso{
  \code{\link{readList}} read in gene set list
}

\examples{
filename=system.file("extdata/gse16873.demo", package = "gage")
demo.data=readExpData(filename, row.names=1)
head(demo.data)
}

\keyword{manip}
