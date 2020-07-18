\name{kegg.gs}
\Rdversion{1.1}
\alias{kegg.gs}
\alias{kegg.gs.dise}
\alias{go.gs}
\alias{carta.gs}
\alias{khier}
\alias{korg}
\alias{bods}
\docType{data}
\title{
  Common gene set data collections
}
\description{
  The gene set data collections derived from KEGG, GO
  and BioCarta databases.
}
\usage{
data(kegg.gs)
data(kegg.gs.dise)
data(go.gs)
data(carta.gs)
}
\format{
  
 kegg.gs is a named list of 177 elements. Each element is a character
 vector of member gene Entrez IDs for a single KEGG pathway. Type
 \code{head(kegg.gs, 3)} for the first 3 gene sets or pathways. Note
 that kegg.gs has been updated since gage version 2.9.1. From then,
 kegg.gs only include the subset of canonical signaling and metabolic
 pathways from KEGG pathway database, and kegg.gs.dise is the subset of
 disease pathways. And it is recommended to do KEGG pathway analysis
 with either kegg.gs or kegg.gs.dise seperately (rather than combined
 altogether) for better defined results. In the near feature, we will
 also generate subsets of go.gs for refined analysis. Note that kegg.gs
 in gageData package still keeps all KEGG pathways, where kegg.gs.sigmet
 and kegg.gs.dise are two subsets of kegg.gs.

 go.gs is a named list of 1000 elements in this package. It is
 a truncated list in this package. The ful list of go.gs has
 ~10000 elements and is provided with an experimental data package
 gageData. Each element is a character
 vector of member gene Entrez IDs for a single Gene Ontology term. Type
 \code{head(go.gs, 3)} for the first 3 gene sets or GO terms.

 carta.gs is a named list of 259 elements. Each element is a character
 vector of member gene Entrez IDs for a single BioCarta pathway. Type
 \code{head(carta.gs, 3)} for the first 3 gene sets or pathways.

 khier is a matrix of 442 rows and 3 columns. This records the category
 and subcategory assignments of all KEGG pathways. The data comes from
 KEGG BRITE Database. This data is used by kegg.gsets function to subset
 the KEGG pathways for more specific pathway analysis.

  korg is a character matrix of ~3000 rows and 6 columns. First 3 columns
  are KEGG species code, scientific name and common name, followed
  columns on gene ID types used for each species: entrez.gnodes ("1" or
  "0", whether EntrezGene is the default gene ID) and representative
  KEGG gene ID and EntrezGene ID. This data comes from pathview package,
  and is used by kegg.gsets function internally.
  
  bods is a character matrix of 19 rows and 4 columns on the mapping
  between gene annotation package names in Bioconductor, common name and
  KEGG code of most common research species. This data comes from
  pathview package, and is used by go.gsets function internally.
}
\details{
These gene set data were compiled using Entrez Gene IDs, gene set names
and mapping information from multiple Bioconductor packages, including:
org.Hs.eg.db, kegg.db, go.db and cMAP. Please check the corresponding
packages for more information.

We only provide gene set data for human with this package. For other
species, please check the experiment data package list of Bioconductor
website or use the Bioconductor package GSEABase to build the users' own
gene set collections.
}
\source{
  Data from multiple Bioconductor packages, including: org.Hs.eg.db,
  kegg.db, go.db and cMAP.
}
\references{
  Entrez Gene <URL: http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=gene>
  KEGG pathways <URL: ftp://ftp.genome.ad.jp/pub/kegg/pathways>
  Gene Ontology <URL: http://www.geneontology.org/>
  cMAP <URL: http://cmap.nci.nih.gov/PW>
  }

\examples{
#load expression and gene set data
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)

data(kegg.gs)
data(go.gs)

#make sure the gene IDs are the same for expression data and gene set
#data
lapply(kegg.gs[1:3],head)
lapply(go.gs[1:3],head)
head(rownames(gse16873))

#GAGE analysis
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis)
gse16873.go.p <- gage(gse16873, gsets = go.gs,
    ref = hn, samp = dcis)
}

\keyword{datasets}
