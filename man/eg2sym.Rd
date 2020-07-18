\name{eg2sym}
\Rdversion{1.1}
\alias{eg2sym}
\alias{sym2eg}
\title{
Conversion between Entrez Gene IDs and official gene
symbols for human genes.
  }
\description{
  These functions convert Entrez Gene IDs to official gene
  symbols for human genes, or vise versa.
}
\usage{
eg2sym(eg)
sym2eg(sym)
}

\arguments{
  \item{eg}{
character vector for Entrez Gene IDs (human genes only).
}
  \item{sym}{
character vector for official gene symbols (human genes only).
}
}
\details{
Currently, only conversion for human genes are supported. Notice that
some gene symbols are not official, hence not recognized and NA will be
returned in such cases.
}
\value{
A character vector giving the converted official gene symbols or Entrez
IDs. 
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
  \code{\link{egSymb}} mapping data between Entrez Gene IDs and official
  symbols;
  \code{\link{readList}} read in gene set list
}

\examples{
#genes in gse16873 was label by Entrez IDs
data(gse16873)
head(rownames(gse16873))
#may convert the gene IDs to official symbols
gse16873.sym<-gse16873
data(egSymb)
rownames(gse16873.sym)<-eg2sym(rownames(gse16873.sym))
head(rownames(gse16873.sym))

#convert kegg.gs correspondingly
data(kegg.gs)
kegg.gs.sym<-lapply(kegg.gs, eg2sym)
lapply(kegg.gs.sym[1:3],head)
#GAGE analysis with the converted data
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
gse16873.kegg.p2 <- gage(gse16873.sym, gsets = kegg.gs.sym,
    ref = hn, samp = dcis)
}

\keyword{manip}
