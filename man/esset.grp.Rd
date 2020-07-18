\name{esset.grp}
\Rdversion{1.1}
\alias{esset.grp}

\title{
The non-redundant signcant gene set list
}
\description{
This function extract a non-redundant signcant gene set list, groups of
redundant gene sets, and related data from \code{gage}
results. Redundant gene sets are those overlap heavily in their
effective member gene lists or core genes.
}
\usage{
esset.grp(setp, exprs, gsets, ref = NULL, samp = NULL, test4up = TRUE,
same.dir = TRUE, compare = "paired", use.fold = TRUE, cutoff = 0.01,
use.q = FALSE, pc = 10^-10, output = TRUE, outname = "esset.grp",
make.plot = FALSE, pdf.size = c(7, 7), core.counts = FALSE, get.essets =
TRUE, bins = 10, bsize = 1, cex = 0.5, layoutType = "circo", name.str =
c(10, 100), ...)
}

\arguments{
    \item{setp}{
    a numeric matrix, the result p-value matrix returned by \code{gage}
    function. Check \code{gage} help information for details.
  }
  \item{exprs}{
    an expression matrix or matrix-like data structure, with genes as
    rows and samples as columns.
  }
  \item{gsets}{
    a named list, each element contains a gene set that is a character
    vector of gene IDs or symbols. For example, type \code{head(kegg.gs)}. A
    gene set can also be a "smc" object defined in PGSEA package.
    Make sure that the same gene ID system is used for both \code{gsets}
    and \code{exprs}.
  }
  \item{ref}{
    a numeric vector of column numbers for the reference condition or
    phenotype (i.e. the control group) in the exprs data matrix. Default
    ref = NULL, all columns are considered as target experiments.
  }
  \item{samp}{
    a numeric vector of column numbers for the target condition or
    phenotype (i.e. the experiment group) in the exprs data
    matrix. Default samp = NULL, all columns other than ref are
    considered as target experiments.
  }
  \item{test4up}{
    boolean, whether the input \code{gage} result or signficant gene sets are
    test results for up-regulated gene sets or not. This information is
    needed for selecting core member genes which contribute to the
    overall signficance of a gene sets.
}
  \item{same.dir}{
    boolean, whether the input \code{gage} result test for changes in a gene set toward a single direction
    (all genes up or down regulated) or changes towards both directions
    simultaneously.
  }
  \item{compare}{
    character, which comparison scheme to be used: 'paired', 'unpaired',
    '1ongroup', 'as.group'. 'paired' is the default, ref and samp are of
    equal length and one-on-one paired by the original experimental
    design; 'as.group', group-on-group comparison between ref and samp;
    'unpaired' (used to be '1on1'), one-on-one comparison between all
    possible ref and samp combinations, although the original
    experimental design may not be one-on-one paired; '1ongroup',
    comparison between one samp column at a time vs the average of all
    ref columns.
  }
  \item{use.fold}{
    Boolean, whether the input \code{gage} results used fold changes or t-test statistics as per
    gene statistics. Default use.fold= TRUE.
  }
  \item{cutoff}{
    numeric, p- or q-value cutoff, between 0 and 1. Default 0.01 (for
    p-value). When q-value is used, recommended cutoff value is 0.1.
}
  \item{use.q}{
    boolean, whether to use q-value or not as the pre-selection of
    a signficant gene set list. Default to be FALSE, i.e. use the
    p-value instead.
}
  \item{pc}{
    numeric, cutoff p-value for the overlap between gene sets to be
    called 'redundant', default to \code{10e-10}, may need
    trial-and-error to find the best value.
  }
  \item{output}{
    boolean, whether output the non-redundant gene set list as
    tab-delimited text file? Default to be TRUE.
  }
  \item{outname}{
    character, the prefix used to label the output file names when
    output = TRUE.
}
  \item{make.plot}{
    boolean, whether to generate the network graph to visualize the
    redundancy (overlap in core genes) between significant gene
    sets. Currently the only feasible option is FALSE.
}
  \item{pdf.size}{
    numeric vector of length 2, spcifies the PDF file size for network
    graph outpout. Currently unsupported.
  }
  \item{core.counts}{
    Currently unsupported.
}
  \item{get.essets}{
    Currently unsupported.
}
  \item{bins}{
    Currently unsupported.
}
  \item{bsize}{
    Currently unsupported.
}
  \item{cex}{
    Currently unsupported.
}
  \item{layoutType}{
    Currently unsupported.
}
\item{name.str}{
  numeric vector of length 2, specifies the substring range of the gene
  set name to show in the network graph. Currently unsupported.
}
  \item{\dots}{
    extra arguments to be passed into internal function
    \code{make.graph}. Currently unsupported.
}
}
\details{
  Redundant gene sets are defined to be those overlap heavily in their effective
  member gene lists or core genes. Core genes are those member genes
  that really contribute to the signficance of the gene set in GAGE
  analysis in the interesting direction(s). Argument \code{pc} set the cutoff for
  the overlap to be called "redundant". The redundancy between
  gene sets is then represented by a undirected graph/network. Groups of
  redundant gene sets are then derived as the connected component in the
  network graph.
  
  The selection criterion for gene sets here is p-value, instead of the
  commonly used q-value. This is because for extracting a non-redundant
  list of signficant gene sets, p-value is relative stable, but q-value
  changes when the total number of gene sets being considered
  changes. Of course, q-value is also a sensible selection criterion,
  when one take this step as a further refinement on the list of
  signficant gene sets.
}
\value{
  The value returned by \code{pairData} is a list of 7 elements:
  \item{essentialSets }{
    character vector, the non-redundant signficant gene set list.
    }
    \item{setGroups }{
      list, each element is a character vector of a group
      of redundant gene sets.
    }
    \item{allSets }{
    character vector, the full list of signficant gene sets.
    }
    \item{setGroups }{
      list, each element is a character vector of a connected component
      in the redundancy graph representation of the gene set.
    }
    \item{overlapCounts }{
      numeric matrix, the overlap core gene counts between the
      signficant gene sets.
    }
    \item{overlapPvals }{
      numeric matrix, the significance (in p-values) of the overlap core gene counts between the
      signficant gene sets.
    }
    \item{coreGeneSets }{
      list, each element is a character vector of the core genes in a
      significant gene set.
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
  \code{\link{gage}} the main function for GAGE analysis;
  \code{\link{sigGeneSet}} significant gene set from GAGE analysis;
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
gse16873.kegg.esg.up <- esset.grp(gse16873.kegg.p$greater,
    gse16873, gsets = kegg.gs, ref = hn, samp = dcis,
    test4up = TRUE, output = TRUE, outname = "gse16873.kegg.up", make.plot = FALSE)
gse16873.kegg.esg.dn <- esset.grp(gse16873.kegg.p$less,
    gse16873, gsets = kegg.gs, ref = hn, samp = dcis,
    test4up = FALSE, output = TRUE, outname = "gse16873.kegg.dn", make.plot = FALSE)
gse16873.kegg.esg.2d <- esset.grp(gse16873.kegg.2d.p$greater,
    gse16873, gsets = kegg.gs, ref = hn, samp = dcis,
    test4up = TRUE, output = TRUE, outname = "gse16873.kegg.2d", make.plot = FALSE)
names(gse16873.kegg.esg.up)
head(gse16873.kegg.esg.up$essentialSets, 4)
head(gse16873.kegg.esg.up$setGroups, 4)
head(gse16873.kegg.esg.up$coreGeneSets, 4)
}

\keyword{htest}
\keyword{multivariate}
\keyword{manip}
