\name{egSymb}
\Rdversion{1.1}
\alias{egSymb}
\docType{data}
\title{
Mapping between Entrez Gene IDs and official symbols
}
\description{
A two column matrix listing the Entrez IDs and official symbols
for all currently known human genes. 
}
\usage{data(egSymb)}
\format{
  The format is:
 chr [1:40784, 1:2] "1" "10" "100" "1000" ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "eg" "sym"
}
\details{
This mapping matrix is commonly used together with functions
\code{eg2sym} and \code{sym2eg}. Check the examples below.
}
\source{
This mapping data matrix were compiled using the gene data from NCBI
Entrez Gene database.

Similar information can also be derived from Bioconductor package
org.Hs.eg.db. Please check the package for more information.
}
\references{
  Entrez Gene <URL: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene>
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
\keyword{datasets}
