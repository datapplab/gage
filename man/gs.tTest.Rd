\name{gs.tTest}
\Rdversion{1.1}
\alias{gs.tTest}
\alias{gs.zTest}
\alias{gs.KSTest}

\title{
  Gene set differential expression test
}
\description{
  These functions test for perturbation of gene sets relative to all genes
  in the microarray data. They are the testing module for \code{gage} and
  single array analysis workflow.

  They use different statistical tests: gs.tTest
  uses two-sample t-test, gs.zTest uses one-sample z-test, gs,KSTest uses
  Kolmogorov-Smirnov test. 
}
\usage{
gs.tTest(exprs, gsets, set.size = c(10, 500), same.dir = TRUE, ...)
gs.zTest(exprs, gsets, set.size = c(10, 500), same.dir = TRUE, ...)
gs.KSTest(exprs, gsets, set.size = c(10, 500), same.dir = TRUE, ...)
}

\arguments{
  \item{exprs}{
    an expression matrix or matrix-like data structure, with genes as
    rows and samples as columns.
  }
  \item{gsets}{
    a named list, each element contains a gene set that is a character
    vector of gene IDs or symbols. For example, type head(kegg.gs). A
    gene set can also be a "smc" object defined in PGSEA package.
    Make sure that the same gene ID system is used for both \code{gsets}
    and \code{exprs}.
  }
  \item{set.size}{
    gene set size (number of genes) range to be considered for
    enrichment test. Tests for too small or too big gene sets are not
    robust statistically or informative biologically. Default to be
    set.size = c(10, 500).
  }
  \item{same.dir}{
    whether to test for changes in a gene set toward a single direction
    (all genes up or down regulated) or changes towards both directions
    simultaneously. For experimentally derived gene sets, GO term
    groups, etc, coregulation is commonly the case, hence same.dir =
    TRUE (default); In KEGG, BioCarta pathways, genes frequently are not
    coregulated, hence it could be informative to let same.dir =
    FALSE. Although same.dir = TRUE could also be interesting for
    pathways.
  }
  \item{\dots}{
    other arguments to be passed into the secondary functions, not used
    currently. 
  }
}
\details{
  These functions are the gene set test module for \code{gage} and
  single array analysis workflow. When used in \code{gage} function, the
  function names are optional values for \code{saaTest} argument. Check
  help information for \code{gage} for details.

  These functions may also used independently without calling \code{gage}
  function.
}

\value{
  As the raw results of gene set tests, a list of 5 components is returned:
  \item{results }{matrix of test statistics, gene sets are rows,
    samp-ref pairs are columns}
  \item{p.results }{matrix of p-values for up-regulation (greater than) tests, gene
    sets are rows, samp-ref pairs are columns}
  \item{ps.results }{matrix of p-values for down-regulation (less than) tests, gene
    sets are rows, samp-ref pairs are columns}
  \item{mstat }{vector of  test statistics mean for individual gene
    sets. Normally, its absoluate value
    measures the magnitude of gene-set level changes, and its sign
    indicates direction of the changes. For gs.KSTest,
    mstat is always positive.}
  \item{setsizes }{vector of effective set size (number of genes) individual gene
    sets}
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
  \code{\link{gage}} the main function for GAGE analysis
}

\examples{
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
data(kegg.gs)

#kegg test
exprs.gage = gagePrep(gse16873, ref = hn, samp = dcis)
str(exprs.gage)
rawRes = gs.tTest(exprs.gage, gsets = kegg.gs)
str(rawRes)
head(rawRes$results)
head(rawRes$p.results)
}

\keyword{htest}
\keyword{multivariate}
