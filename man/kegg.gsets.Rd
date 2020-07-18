\name{kegg.gsets}
\alias{kegg.gsets}
\title{
  Generate up-to-date KEGG pathway gene sets
}
\description{
  Generate up-to-date KEGG pathway gene sets for any specified
  KEGG species.
}
\usage{
kegg.gsets(species = "hsa", id.type = "kegg", check.new=FALSE)
}
\arguments{
  \item{species}{
    character, either the kegg code, scientific name or the common name of
    the target species. This applies to both pathway and gene.data or
    cpd.data. When KEGG ortholog pathway is considered,
    species="ko". Default species="hsa", it is equivalent to use either
    "Homo sapiens" (scientific name) or "human" (common name). 
  }
  \item{id.type}{
    character, desired ID type for the get sets, case insensitive. Default
    idtype="kegg", i.e. the primary KEGG gene ID. The other valid
    option is "entrez", i.e. Entrez Gene. Entrez Gene is the primary KEGG gene ID
    for many common model organisms, like human, mouse, rat etc, hence these
    two options have the same effect. For other species, primary KEGG gene
    ID is not  Entrez Gene. 
  }
  \item{check.new}{
    logical, whether to check for any newly added reference pathways
    online. Default to FALSE. Because such online checking takes time
    and new reference pathways are rarely added, it is usually not
    suggested to set this argument TRUE, unless this is really needed.
  }
}
\details{
  The latest KEGG pathway gene sets are derived by connecting to the
  database in real time. This way, we can create high quality gene set
  data for pathway analysis for over 2400 KEGG species.
  
  Note that we have generated GO gene set for 4 species, human, mouse,
  rat, yeast as well as KEGG Ortholog, and provided the data in package
  gageData.
}
\value{
  A named list with the following elements:
  \item{kg.sets}{
    KEGG gene sets, a named list. Each element is a character
    vector of member gene IDs for a single KEGG pathway. The number of
    elements of this list is the total number of KEGG pathways defined for
    the specified species.
  }
  \item{sigmet.idx}{
    integer indice, which elements in kg.sets are signaling or metabolism pathways.
  }

  \item{sig.idx}{
    integer indice, which elements in kg.sets are signaling pathways.
  }

  \item{met.idx}{
    integer indice, which elements in kg.sets are metabolism pathways.
  }

  \item{dise.idx}{
    integer indice, which elements in kg.sets are disease pathways.
  }
  The *.idx elements here are all used to subset kg.sets for more specific
  type pathway aanlysis.
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
  \code{\link{kegg.gs}} for precompiled KEGG and other common gene set
  data collection
}
\examples{
#GAGE analysis use the latest KEGG pathway  definitions, instead of
#kegg.gs
kg.hsa=kegg.gsets()
data(gse16873)
hn=(1:6)*2-1
dcis=(1:6)*2
kegg.sigmet=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
gse16873.kegg.p <- gage(gse16873, gsets = kegg.sigmet,
                        ref = hn, samp = dcis)

#E coli KEGG Id is different from Entre Gene
kg.eco=kegg.gsets("eco")
kg.eco.eg=kegg.gsets("eco", id.type="entrez")
head(kg.eco$kg.sets,2)
head(kg.eco.eg$kg.sets,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
