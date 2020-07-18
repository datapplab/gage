\name{sigGeneSet}
\Rdversion{1.1}
\alias{sigGeneSet}

\title{
Significant gene set from GAGE analysis
}
\description{
  This function sorts and counts signcant gene sets based on q- or p-value
  cutoff.
}
\usage{
sigGeneSet(setp, cutoff = 0.1, dualSig = (0:2)[2], qpval = c("q.val",
"p.val")[1],heatmap=TRUE, outname="array", pdf.size = c(7,7),
p.limit=c(0.5, 5.5), stat.limit=5,  ...)
}

\arguments{
  \item{setp}{
    the result object returned by \code{gage} function, either a numeric
    matrix or a list of two such matrices. Check \code{gage} help
    information for details.
}
  \item{cutoff}{
    numeric, q- or p-value cutoff, between 0 and 1. Default 0.1 (for
    q-value). When p-value is used, recommended cutoff value is 0.001 for data with more than 2
    replicates per condition or 0.01 for les sample sizes.
}
  \item{dualSig}{
    integer, switch argument controlling how dual-signficant gene sets
    should be treated. This argument is only useful when Stouffer method
    is not used in gage function (use.stouffer=FALSE), hence makes no
    difference normally. 0: discard such gene sets from
    the final significant gene set list; 1: keep such gene sets in the
    more significant direction and remove them from the less significant
    direction; 2: keep such gene sets in the lists for both
    directions. default to 1. Dual-signficant means a gene set is called
    significant simultaneously in both 1-direction tests (up- and
    down-regulated). Check the details for more information.
}
  \item{qpval}{
    character, specifies the column name used for gene set selection,
    i.e. what type of q- or p-value to use in gene set
    selection. Default to 
    be "q.val" (q-value using BH procedure). "p.val" is the unadjusted
    global p-value and may be used as selection criterion
    sometimes. 
  }
  \item{heatmap}{
    boolean, whether to plot heatmap for the selected gene data as a PDF
    file. Default to be FALSE.
  }
  \item{outname}{
    a character string, to be used as the prefix of the output data
    files. Default to be "array".
  }
  \item{pdf.size}{
    a numeric vector to specify the  the width and height of PDF
    graphics region in inches. Default to be c(7, 7).
  }
  \item{stat.limit}{
    numeric vector of length 1 or 2 to specify the value range of gene
    set statistics to visualize using the heatmap. Statistics beyong will be
    reset to equal the proximal limit. Default to 5, i.e. plot all gene set statistics
    within (-5, 5) range. May also be NULL, i.e. plot all statistics
    without limit. This argument allows optimal differentiation between most gene
    set statistic values when extremely positive/negative values exsit and squeeze
    the normal-value region.
  }
  \item{p.limit}{
    numeric vector of length 1 or 2 to specify the value range of gene
    set -log10(p-values) to visualize using the heatmap. Values beyong will be
    reset to equal the proximal limit. Default to c(0.5,5.5), i.e. plot all -log10(p-values)
    within this range. This argument is similar to argument stat.limit.
  }
  \item{\dots}{
    other arguments to be passed into the inside \code{gs.heatmap}
    function, which is a wrapper of the \code{heatmap2} function.
  }
}

\details{
By default, heatmaps are produced to show the gene set
perturbations using either -log10(p-value) or statistics.

Since gage package version 2.2.0, Stouffer's method is used as the
  default procedure for more robust p-value summarization. With the original p-value
  summarization, i.e. negative log sum following a Gamma distribution as
  the Null hypothesis, the global p-value could be
  heavily affected by a small subset of extremely small individual
  p-values from pair-wise comparisons. Such sensitive global p-value
  leads to the "dual signficance" phenomenon. In other words, Gene sets
  are signficantly up-regulated in a subset 
    of experiments, but down-regulated in another subset. Note that dual-signficant gene sets are not the
    same as gene sets called signficant in 2-directional tests, although
    they are related.
}
\value{
  \code{sigGeneSet} function returns a named list of the same structure
  as \code{gage} result. Check \code{gage} help information for
  details. 
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
  \code{\link{esset.grp}} non-redundant signcant gene set list;
  \code{\link{essGene}} essential member genes in a gene set;
}

\examples{
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
data(kegg.gs)

#kegg test for 1-directional changes
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis)
#kegg test for 2-directional changes
gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis, same.dir = FALSE)
gse16873.kegg.sig<-sigGeneSet(gse16873.kegg.p, outname="gse16873.kegg")
str(gse16873.kegg.sig)
gse16873.kegg.2d.sig<-sigGeneSet(gse16873.kegg.2d.p, outname="gse16873.kegg")
str(gse16873.kegg.2d.sig)
#also check the heatmaps in pdf files named "*.heatmap.pdf".
}

\keyword{htest}
\keyword{multivariate}
\keyword{manip}
