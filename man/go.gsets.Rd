\name{go.gsets}
\alias{go.gsets}
\title{
  Generate up-to-date GO (Gene Ontology) gene sets
}
\description{
  Generate up-to-date GO gene sets for specified species, which has
  either Bioconductor or user supplied gene annotation package.
}
\usage{
go.gsets(species = "human", pkg.name=NULL, id.type = "eg", keep.evidence=FALSE)
}
\arguments{
  \item{species}{
    character, common name of the target species. Note that scientific
    name is not used here. Default species="human". For a list of other
    species names, type in: {data(bods); bods}
  }
  \item{pkg.name}{
    character, the gene annotation package name for the target
    species. It is either one of the Bioconductor OrgDb packages or a
    user supplied package with similar data structures created using
    AnnotationForge package. Default species=NULL, the right annotation
    package will be located for the specified speices in Bioconductor
    automatically. Otherwise, the user need to prepare and supply a
    usable annotation package in the same format.
  }
  \item{id.type}{
    character, target ID type for the get sets, case insensitive. Default
    gene.idtype="eg", i.e. Entrez Gene. Entrez Gene is the primary GO gene ID
    for many common model organisms, like human, mouse, rat etc. For
    other species may use "orf" (open reading frame) and other specific
    gene IDs, please type in: {data(bods); bods} to check the details. 
  }
  \item{keep.evidence}{
    boolean, whether to keep the GO evidence codes for genes assigned to
    each GO term. The evidence codes describe what type of evidence was
    present in that reference to make the annotation. Default
    keep.evidence = FALSE, which removes all the evidence codes and
    redundant gene entries due to multiple evidence codes in a gene set.
  }
}
\details{
  The updated GO gene sets are derived using Bioconductor or user
  supplied gene annotation packages. This way, we can create gene set
  data for GO analysis for 19 common species annotated in Bioconductor
  and more others by the users. Note that we have generated GO gene set
  for 4 species, human, mouse, rat and yeast, and provided the data in
  package gageData.
}
\value{
  A named list with the following elements:
  \item{go.sets}{
    GO gene sets, a named list. Each element is a character
    vector of member gene IDs for a single GO. The number of
    elements of this list is the total number of GO terms defined for
    the specified species.
  }
  \item{go.subs}{
    go.subs is a named list of three elements: "BP", "CC" and "MF",
    corresponding to biological process, cellular component and molecular
    function subtrees. It may be more desirable to conduct separated GO
    enrichment test on each of these 3 subtrees as shown in the example code.
  }
  \item{go.mains}{
    go.mains is a named vector of three elements: "BP", "CC" and "MF",
    corresponding to the root node of biological process, cellular component and molecular
    function subtree, i.e. their corresponding indecies in go.sets.
  }
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
  \code{\link{go.gs}} for precompiled GO and other common gene set
  data collection
}
\examples{
#GAGE analysis use the updated GO  definitions, instead of
#go.gs.
#the following lines are not run due to long excution time.
#The data generation steop may take 10-30 seconds. You may want
#save the gene set data, i.e. go.hs, for future uses.
#go.hs=go.gsets(species="human")
#data(gse16873)
#hn=(1:6)*2-1
#dcis=(1:6)*2
#go.bp=go.hs$go.sets[go.hs$go.subs$BP]
#gse16873.bp.p <- gage(gse16873, gsets = go.bp,
#                        ref = hn, samp = dcis)

#Yeast and few othr species gene Id is different from Entre Gene
data(bods)
bods
#not run
#go.sc=go.gsets("Yeast")
#lapply(go.sc$go.sets[1:3], head, 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
