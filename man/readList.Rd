\name{readList}
\Rdversion{1.1}
\alias{readList}

\title{
Read in gene set data as a named list
}
\description{
  This function reads in gene set data in GMT (.gmt) format as a named
  list. GMT is defined originally by GSEA program. The code may be
  slightly revised for reading in gene set data in other tab-delimited formats
  too.
}
\usage{
readList(file)
}

\arguments{
  \item{file}{
character string, the full path name to the gene set data file in
GMT format.
}}

\value{
A named list, each element is a character vector giving the gene IDs of a
gene set.
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
  \code{\link{readExpData}} read in gene expression data
}

\examples{
#an example GMT gene set data derived from MSigDB data
filename=system.file("extdata/c2.demo.gmt", package = "gage")
demo.gs=readList(filename)
demo.gs[1:3]
#to use these gene sets with gse16873, need to convert the gene symbols
#to Entrez IDs first
data(egSymb)
demo.gs.sym<-lapply(demo.gs, sym2eg)
demo.gs.sym[1:3]
}

\keyword{manip}
