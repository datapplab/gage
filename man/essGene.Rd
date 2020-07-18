\name{essGene}
\Rdversion{1.1}
\alias{essGene}

\title{
Essential member genes in a gene set
}
\description{
This function extracts data for essential member genes in a gene
set. Essential genes are genes that have changes over noise level.
}
\usage{
essGene(gs, exprs, ref = NULL, samp = NULL, gsets = NULL, compare
= "paired", use.fold = TRUE, rank.abs = FALSE, use.chi = FALSE, chi.p =
0.05, ...)
}

\arguments{
  \item{gs}{
character, either the name of an interesting gene set in a gene set
collection passed by \code{gsets} argument, or a vector of gene IDs.
    Make sure that the same gene ID system is used for both \code{gs}
    and \code{exprs}.
}
  \item{exprs}{
    an expression matrix or matrix-like data structure, with genes as
    rows and samples as columns.
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
  \item{gsets}{
    a named list, each element contains a gene set that is a character
    vector of gene IDs or symbols. For example, type \code{head(kegg.gs)}. A
    gene set can also be a "smc" object defined in PGSEA package.
    Make sure that the same gene ID system is used for both \code{gsets}
    and \code{exprs}. Default to be NULL, then argument \code{gs} needs to be a
    vector of gene IDs.
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
  \item{rank.abs}{
boolean, whether to sort the essential gene data based on absoluate
changes. Default to be FALSE.
}
  \item{use.chi}{
boolean, whether to use chi-square test to select the essential genes.
Default to be FALSE, use the mean plus standard deviation of all gene
changes instead. Check details for more information.
}
  \item{chi.p}{
numeric value between 0 and 1, cutoff p-value for the chi-square test to
select the essential genes. Default to 0.05.
}  \item{\dots}{
    other arguments to be passed into the inside \code{gagePrep} function.
  }
}
\details{
There are two different criteria for essential gene selection. One uses
a chi-square test to determin whether the change of a gene is more than
noise. A second considers any changes beyond 1 standard deviation from
mean of all genes as real.

Note that essential genes are different from core genes considered in
\code{esset.grp} function. Essential genes may change in a different
direction than the overall change of a gene set. But core genes need to
change in the in the interesting direction(s) of the gene set test.
}
\value{
  A expression data matrix extracted for the essential genes, with
  similar structure as \code{exprs}.
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
  \code{\link{geneData}} output and visualization of expression data
  for selected genes;
  \code{\link{esset.grp}} non-redundant signcant gene set list;
}

\examples{
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)

#kegg test for 1-directional changes
data(kegg.gs)
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis)
rownames(gse16873.kegg.p$greater)[1:3]
gs=unique(unlist(kegg.gs[rownames(gse16873.kegg.p$greater)[1:3]]))
essData=essGene(gs, gse16873, ref =hn, samp =dcis)
head(essData)
ref1=1:6
samp1=7:12
#generated text file for data table, pdf files for heatmap and scatterplot
for (gs in rownames(gse16873.kegg.p$greater)[1:3]) {
    outname = gsub(" |:|/", "_", substr(gs, 10, 100))
    geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
        samp = samp1, outname = outname, txt = TRUE, heatmap = TRUE,
        Colv = FALSE, Rowv = FALSE, dendrogram = "none", limit = 3, scatterplot = TRUE)
}
}

\keyword{multivariate}
\keyword{manip}
