\name{heter.gage}
\Rdversion{1.1}
\alias{heter.gage}
\alias{pairData}

\title{
  GAGE analysis for heterogeneous data
}
\description{
  \code{heter.gage} is a wrapper function of \code{gage} for heterogeneous
  data. \code{pairData} prepares the heterogeneous data and related
  arguments for GAGE analysis.
}
\usage{
heter.gage(exprs, gsets, ref.list, samp.list, comp.list = "paired",
use.fold = TRUE, ...)

pairData(exprs, ref.list, samp.list, comp.list = "paired", use.fold =
TRUE, ...)
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
  \item{ref.list}{
    a list of \code{ref} inputs for \code{gage} function. In other
    words, each element of the list is a column number vector for the reference condition or
    phenotype (i.e. the control group) in the exprs data matrix.
  }
  \item{samp.list}{
    a list of \code{samp} inputs for \code{gage} function. In other
    words, each element of the list is a column number vector for the target condition or
    phenotype (i.e. the experiment group) in the exprs data
    matrix.
  }
  \item{comp.list}{
    a list or a vector of \code{compare} input(s) for \code{gage}
    function. The length of the list or vector should equal to the length
    of \code{ref.list} and \code{samp.list} or 1. In the latter case,
    all analyses will use the same comparison scheme. The same as 
    \code{compare}, the element value(s) in \code{comp.list} can be
    'paired', 'unpaired', '1ongroup' or 'as.group'. Default to be 'paired'.
  }
  \item{use.fold}{
    Boolean, whether to use fold changes or t-test statistics as per
    gene statistics. Default use.fold= TRUE.
  }
  \item{\dots}{
    other arguments to be passed into \code{gage}.
  }
}
\details{
  comp.list can be a list or vector of mixture values of 'paired' and
  'unpaired' matching the experiment layouts of the heterogeneous data. In
  such cases, each ref-samp pairs and corresponding columns in the result
  data matrix after calling \code{pairData} are assigned different weights
  when calling \code{gage} in the next step.
  The inclusion of '1ongroup' and 'as.group' in comp.list would make
  weight assignment very complicated especially when the sample sizes are
  different for the individual experiments of the heterogeneous data.
}
\value{
  The output of \code{pairData} is a list of 2 elements:
  \item{exprs }{a data matrix derived from the input expression data matrix
    \code{exprs}, but ready for column-wise gene est tests. In the
    matrix, genes are rows, and columns are the per gene test
    statistics from the ref-samp pairwise comparison.}
  \item{weights }{weights assigned to columns of the output data matrix
    \code{exprs} when calling \code{gage} next. The value may be NULL if
    \code{comp.list} are all 'paired'.}

  The result returned by \code{heter.gage} function is the same as
  result of \code{gage}, i.e. either a single data matrix
  (same.dir = FALSE, test for two-directional changes) or
  a named list of two data matrix (same.dir = TRUE, test for single-direction
  changes) for the results of up- ($greater) and down- ($less) regulated
  gene sets. Check help information for \code{gage} for details.
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
  \code{\link{gage}} the main function for GAGE analysis;
  \code{\link{gagePipe}} pipeline for multiple GAGE analysis in a batch
}

\examples{
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
data(kegg.gs)

library(gageData)
data(gse16873.2)
cn2=colnames(gse16873.2)
hn2=grep('HN',cn2, ignore.case =TRUE)
dcis2=grep('DCIS',cn2, ignore.case =TRUE)

#combined the two half dataset
gse16873=cbind(gse16873, gse16873.2)
refList=list(hn, hn2+12)
sampList=list(dcis, dcis2+12)

#quick look at the heterogeneity of the combined data
summary(gse16873[,hn[c(1:2,7:8)]])
#if graphic devices open:
#boxplot(data.frame(gse16873))
gse16873.kegg.heter.p <- heter.gage(gse16873, gsets = kegg.gs,
    ref.list = refList, samp.list = sampList)
gse16873.kegg.heter.2d.p <- heter.gage(gse16873, gsets = kegg.gs,
    ref.list = refList, samp.list = sampList, same.dir = FALSE)
str(gse16873.kegg.heter.p)
head(gse16873.kegg.heter.p$greater[, 1:5])
}

\keyword{htest}
\keyword{multivariate}
\keyword{manip}
