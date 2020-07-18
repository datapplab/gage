\name{gagePipe}
\Rdversion{1.1}
\alias{gagePipe}

\title{
GAGE analysis pipeline
}
\description{
  Function gagePipe runs mutliple rounds of GAGE in a
  batch without interference, and outputs signcant gene set lists in text
  format, heatmaps in pdf format, and save the results in RData format.
}
\usage{
gagePipe(arraydata, dataname = "arraydata", trim.at=TRUE, sampnames, gsdata = NULL,
gsname = c("kegg.gs", "go.gs"), ref.list, samp.list, weight.list = NULL,
comp.list = "paired", q.cutoff = 0.1, heatmap=TRUE, pdf.size = c(7,
7), p.limit=c(0.5, 5.5), stat.limit=5, ... )
}

\arguments{
  \item{arraydata}{
    corresponds to \code{exprs} argument for \code{gage} function. But
    can either be a matrix-like data structure when the data has been
    loaded into R or the full path to the data file in .RData format if
    the data has not been loaded.
  }
  \item{dataname}{
    character, the name of the data to be analyzed. This name will be
    included as the prefix of the output file names. Default to be
    "arraydata".
  }
  \item{trim.at}{
    boolean, whether to trim the suffix "_at" from the probe set IDs or
    row names of the microarray data. Default to be TRUE.
  }
  \item{sampnames}{
    character vector, the names of the sample groups, on which the GAGE
    analysis is done. Each sample groups corresponds to one element of
    samp.list and the matching element of ref.list. These names will be 
    included in the output file names or object names.
  }
  \item{gsdata}{
    character, the full path to the gene set data file in .RData format if
    the data has not been loaded. Default to be NULL, i.e. the gene set
    data has been loaded.
    Make sure that the same gene ID system is used for both \code{gsdata}
    and \code{arraydata}.
}
  \item{gsname}{
    character, the name(s) of the gene set collections to be
    used. Default to be \code{c("kegg.gs", "go.gs")}.
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
  \item{weight.list}{
    a list or a vector of \code{weights} input(s) for \code{gage}
    function. As a list, the length of \code{weight.list} should equal to the length
    of \code{ref.list} and \code{samp.list} or 1.  The
    \code{weight.list} vector or its element
    vectors of  should match the elements of \code{ref.list} and
    \code{samp.list} in length or being NULL. When \code{weight.list} is
    a vector or length 1 list, 
    all analyses will use the same \code{weights} setting.
  }
  \item{comp.list}{
    a list or a vector of \code{compare} input(s) for \code{gage}
    function. The length of the list or vector should equal to the length
    of \code{ref.list} and \code{samp.list} or 1. In the latter case,
    all analyses will use the same comparison scheme. The same as 
    \code{compare}, the element value(s) in \code{comp.list} can be
    'paired', 'unpaired', '1ongroup' or 'as.group'. Default to be 'paired'.
  }
  \item{q.cutoff}{
    numeric, q-value cutoff between 0 and 1 for signficant gene sets
    selection. Default to be 0.1.
  }
  \item{heatmap}{
    boolean, whether to plot heatmap for the selected gene data as a PDF
    file. Default to be FALSE.
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
    other arguments to be passed into \code{gage} or \code{gs.heatmap}
    function, which is a wrapper of the \code{heatmap2} function.
  }
}
\details{
  \code{gagePipe} carries two rounds of GAGE analysis on each sample
  groups for each
  gene set collection specified in \code{gsnames}: one test for
  1-direction changes (up- or down-regualted gene sets), one test for
  2-direction changes (two-way perturbed gene sets). Correspondingly,
  the \code{gage} result p-value matrices for the signficant gene sets are written
  into two tab-delimited text files, named after the \code{dataname}
  and \code{sampnames}. Note that the text file for 1-direction changes
  tests combines results for both up- and down-regulated gene
  sets. By default, heatmaps in pdf format are also produced to show the gene set
  perturbations using either -log10(p-value) or statistics. Meanwhile,
  the full \code{gage} analysis result objects 
  (named lists of p-value or statistics matrices) are saved into a .RData file. The
  result objects are name after the \code{sampnames} and
  \code{gsnames}.
}
\value{
  The function returns invisible 1 when successfully
  executed. 
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
  \code{\link{heter.gage}} GAGE analysis for heterogeneous data
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

#multiple GAGE analysis in a batch with the combined data
gse16873=cbind(gse16873, gse16873.2)
dataname='gse16873' #output data prefix
sampnames=c('dcis.1', 'dcis.2')
refList=list(hn, hn2+12)
sampList=list(dcis, dcis2+12)
gagePipe(gse16873, gsname = "kegg.gs", dataname = "gse16873",
    sampnames = sampnames, ref.list = refList, samp.list = sampList,
    comp.list = "paired")

#follow up comparison between the analyses
load('gse16873.gage.RData')
#list gage result objects
objects(pat = "[.]p$")
gageComp(sampnames, dataname, gsname = "kegg.gs",
    do.plot = TRUE)
}

\keyword{htest}
\keyword{multivariate}

