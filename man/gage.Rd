\name{gage}
\Rdversion{1.1}
\alias{gage}
\alias{gagePrep}
\alias{gageSum}

\title{
  GAGE (Generally Applicable Gene-set Enrichment) analysis
}
\description{
  Run GAGE analysis to infer gene sets (or pathways, functional groups
  etc) that are signficantly perturbed relative to all genes
  considered. GAGE is generally applicable to essentially all microarray
  dta independent of data attributes including sample size, experimental
  layout, study design, and all types of heterogeneity in data
  generation.
  
  gage is the main function; gagePrep is the functions for the initial
  data preparation, especially sample pairing; gageSum carries out the
  final meta-test summarization.
}

\usage{
gage(exprs, gsets, ref = NULL, samp = NULL, set.size = c(10, 500),
same.dir = TRUE, compare = "paired", rank.test = FALSE, use.fold = TRUE,
FDR.adj = TRUE, weights = NULL, full.table = FALSE, saaPrep = gagePrep,
saaTest = gs.tTest, saaSum = gageSum, use.stouffer=TRUE, ...)

gagePrep(exprs, ref = NULL, samp = NULL, same.dir = TRUE, compare =
"paired", rank.test = FALSE, use.fold = TRUE, weights = NULL, full.table =
FALSE, ...)

gageSum(rawRes, ref = NULL, test4up = TRUE, same.dir =
TRUE, compare = "paired", use.fold = TRUE, weights = NULL, full.table =
FALSE, use.stouffer=TRUE, ...)
}

\arguments{
  \item{exprs}{
    an expression matrix or matrix-like data structure, with genes as
    rows and samples as columns.
  }
  \item{gsets}{
    a named list, each element contains a gene set that is a character
    vector of gene IDs or symbols. For example, type \code{head(kegg.gs)}. A
    gene set can also be a "smc" object defined in PGSEA package. Please
    make sure that the same gene ID system is used for both \code{gsets}
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
  \item{set.size}{
    gene set size (number of genes) range to be considered for
    enrichment test. Tests for too small or too big gene sets are not
    robust statistically or informative biologically. Default to be
    set.size = c(10, 500).
  }
  \item{same.dir}{
    boolean, whether to test for changes in a gene set toward a single direction
    (all genes up or down regulated) or changes towards both directions
    simultaneously. For experimentally derived gene sets, GO term
    groups, etc, coregulation is commonly the case, hence same.dir =
    TRUE (default); In KEGG, BioCarta pathways, genes frequently are not
    coregulated, hence it could be informative to let same.dir =
    FALSE. Although same.dir = TRUE could also be interesting for
    pathways.
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

    For PAGE-like analysis, the default is compare='as.group', which is
    the only option provided in the original PAGE method. All other
    comparison schemas are set here for direct comparison to gage.
  }
  \item{rank.test}{
    rank.test: Boolean, whether do the optional rank based two-sample
    t-test (equivalent to the non-parametric Wilcoxon Mann-Whitney test)
    instead of parametric two-sample t-test. Default rank.test =
    FALSE. This argument should be used with respect to argument
    saaTest.
  }
  \item{use.fold}{
    Boolean, whether to use fold changes or t-test statistics as per
    gene statistics. Default use.fold= TRUE.
  }
  \item{FDR.adj}{
    Boolean, whether to do adjust for multiple testing as to control FDR
    (False dicovery rate). Default to be TRUE.
  }
  \item{weights}{
    a numeric vector to specify the weights assigned to pairs of
    ref-samp. This is needed for data with both technical replicates and
    biological replicates as to count for the different contributions
    from the two types of replicates. This argument is also useful in
    manually paring ref-samp for unpaired data, as in pairData function.
    function. Default to be NULL.
  }
  \item{full.table}{
    This option is obsolete. Boolean, whether to output the full table of all individual p-values
    from the pairwise comparisons of ref and samp. Default to
    be FALSE.
  }
  \item{saaPrep}{
    function used for data preparation for single array based analysis,
    including sanity check, sample pairing, per gene statistics
    calculation etc. Default to be gagePrep, i.e. the default data
    preparation routine for gage analysis.
  }
  \item{saaTest}{
    function used for gene set tests for single array based
    analysis. Default to be gs.tTest, which features a two-sample t-test
    for differential expression of gene sets. Other options includes:
    gs.zTest, using one-sample z-test as in PAGE, or gs.KSTest, using
    the non-parametric Kolmogorov-Smirnov tests as in GSEA. The two
    non-default options should only be used when rank.test = FALSE.
  }
  \item{saaSum}{
    function used for summarization of the results from single array
    analysis (i.e. pairwise comparison between ref and samp). This
    function should include a meta-test for a global p-value or summary
    statistis and a FDR adjustment for multi-testing issue. Default to
    be gageSum, i.e. the default data summarization routine for gage
    analysis.
  }
  \item{rawRes}{
    a named list, the raw results of gene set tests. Check the help
    information of gene set test functions like \code{gs.tTest} for
    details. 
}
    \item{test4up}{
    boolean, whether to summarize the p-value results for up-regulation
    test (p.results) or not (ps.results for down-regulation).
    This argument is only needed when
    the argument same.dir=TRUE in the main \code{gage} function,
    i.e. when test for one-directional changes.
}

  \item{use.stouffer}{
    Boolean, whether to use Stouffer's method when summarizing
    individual p-values into a global p-value. Default to
    be TRUE. This default setting is recommended as to avoid the "dual
    significance", i.e. a gene set is significant for both up-regulation
    and down-regulation tests simultaneously. Dual signficance  occurs
    sometimes for data with large sample size, due to extremely small
    p-values in a few pair-wise comparison. More robust p-value
    summarization using Stouffer's method is a
    important new feature added to GAGE since version 2.2.0
    (Bioconductor 2.8). This new argument is set as to provide a option
    to the original summarization based on Gamma distribution (FALSE).
  }
  \item{\dots}{
    other arguments to be passed into the optional functions for
    saaPrep, saaTest and saaSum.
  }
}
\details{
  We proposed a single array analysis (i.e. the one-on-one comparison)
  approach with GAGE. Here we made single array analysis a general
  workflow for gene set analysis.  Single array analysis has 4 major
  steps: Step 1 sample pairing, Step 2 per gene tests, Step 3 gene set
  tests and Step 4 meta-test summarization. Correspondingly, this new main
  function, gage, is divided into 3 relatively independent modules. Module
  1 input preparation covers step 1-2 of single array analysis. Module 2
  corresponds to step 3 gene set test, and module 3 to step 4 meta-test
  summarization. These 3 modules become 3 argument functions to gage,
  saaPrep, saaTest and saaSum. The modulization made gage open to
  customization or plug-in routines at each steps and fully realize the
  general applicability of single array analysis. More examples will be
  included in a second vignette to demonstrate the customization with
  these modules.

  some important updates has been made to gage package since version 2.2.0
  (Bioconductor 2.8):
  First,  more robust p-value summarization using Stouffer's method
  through argument use.stouffer=TRUE. The original p-value
  summarization, i.e. negative log sum following a Gamma distribution as
  the Null hypothesis, may produce less stable global p-values for large
  or heterogenous datasets. In other words, the global p-value could be
  heavily affected by a small subset of extremely small individual
  p-values from pair-wise comparisons. Such sensitive global p-value
  leads to the "dual signficance" phenomenon. Dual-signficant means a gene set is called
    significant simultaneously in both 1-direction tests (up- and
    down-regulated). "Dual
    signficance" could be informative in revealing the sub-types or
    sub-classes in big clinical or disease studies, but may not be
    desirable in other cases.
    Second, output of gage function now includes the gene set test
    statistics from pair-wise comparisons for all proper gene sets. The
    output is always a named list now, with either 3 elements
    ("greater", "less", "stats") for one-directional test or 2 elements
    ("greater", "stats") for two-directional test. 
    Third, the individual p-value (and test statistics)from dependent pair-wise
    comparisions, i.e. comparisions between the same experiment vs
    different controls, are now summarized into a single value. In other
    words, the column number of individual p-values or statistics is
    always the same as the sample number in the experiment (or disease)
    group. This change made the argument value compare="1ongroup"
    and argument full.table less useful. It also became easier to check the
    perturbations at gene-set level for individual samples. 
    Fourth, whole gene-set level changes (either p-values or statistics)
    can now be visualized using heatmaps due to the third change above.
    Correspondingly, functions \code{sigGeneSet} and \code{gagePipe} have
    been revised to plot heatmaps for whole gene sets. 
  }
\value{
  The result returned by gage function is a named list, with either 3 elements
    ("greater", "less", "stats") for one-directional test (same.dir =
    TRUE) or 2 elements ("greater", "stats") for two-directional test
    (same.dir = FALSE). Elements "greater" and "less" are two data
    matrices of the same structure, mainly the p-values, element "stats"
    contains the test statistics. Each data matrix here has gene sets as rows sorted by
    global p- or q-values. Test signficance or statistics columns include:
  \item{p.geomean }{geometric mean of the individual p-values from
    multiple single array based gene set tests}
  \item{stat.mean }{mean of the individual statistics from multiple
    single array based gene set tests. Normally, its absoluate value
    measures the magnitude of gene-set level changes, and its sign
    indicates direction of the changes. When saaTest=gs.KSTest,
    stat.mean is always positive.}
  \item{p.val }{gloal p-value or summary of the individual p-values from
    multiple single array based gene set tests. This is the default
    p-value being used.}
  \item{q.val }{FDR q-value adjustment of the global p-value using
    the Benjamini & Hochberg procedure implemented in multtest
    package. This is the default q-value being used.}
  \item{set.size }{the effective gene set size, i.e. the number of genes
    included in the gene set test}
  \item{other columns }{columns of the individual p-values or
    statistics, each measures the gene set perturbation in a single
    experiment (vs its control or all controls, depends on the
    "compare argument value)}

  The result returned by \code{gagePrep} is a data matrix derived from 
  \code{exprs}, but ready for column-wise gene est tests. In the
  matrix, genes are rows, and columns are the per gene test
  statistics from the ref-samp pairwise comparison.

  The result returned by \code{gageSum} is almost identical to the results
  of \code{gage} function, it is also a named list but has only 2
  elements, "p.glob" and "results", with one round of test results. 
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
  \code{\link{gs.tTest}}, \code{\link{gs.zTest}}, and
  \code{\link{gs.KSTest}} functions used for gene set test;
  \code{\link{gagePipe}} and \code{\link{heter.gage}} function used for
  multiple GAGE analysis in a batch or combined GAGE analysis on
  heterogeneous data
}

\examples{
data(gse16873)
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
data(kegg.gs)
data(go.gs)

#kegg test for 1-directional changes
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis)
#go.gs with the first 1000 entries as a fast example.
gse16873.go.p <- gage(gse16873, gsets = go.gs, 
    ref = hn, samp = dcis)
str(gse16873.kegg.p)
head(gse16873.kegg.p$greater)
head(gse16873.kegg.p$less)
head(gse16873.kegg.p$stats)
#kegg test for 2-directional changes
gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis, same.dir = FALSE)
head(gse16873.kegg.2d.p$greater)
head(gse16873.kegg.2d.p$stats)

###alternative ways to do GAGE analysis###
#with unpaired samples
gse16873.kegg.unpaired.p <- gage(gse16873, gsets = kegg.gs, 
    ref = hn, samp = dcis, compare = "unpaired")

#other options to tweak includes:
#saaTest, use.fold, rank.test, etc. Check arguments section above for
#details and the vignette for more examples.
}

\keyword{htest}
\keyword{multivariate}
\keyword{manip}
