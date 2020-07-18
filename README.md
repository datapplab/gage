
# gage package

[![](https://img.shields.io/badge/release%20version-3.38.3-blue.svg)](https://www.bioconductor.org/packages/gage)
[![](https://img.shields.io/badge/devel%20version-3.39.3-green.svg)](https://github.com/datapplab/gage)
[![](https://img.shields.io/badge/BioC%20since-2009-blue.svg)](https://www.bioconductor.org/packages/devel/bioc/html/gage.html#since)
[![](https://img.shields.io/badge/GitHub%20since-2020-green.svg)](https://github.com/datapplab/gage)

## Overview

GAGE is a widely used method for gene set (enrichment or GSEA) or pathway analysis. GAGE is generally applicable independent of microarray or RNA-Seq data attributes including sample sizes, experimental designs, assay platforms, and other types of heterogeneity, and consistently achieves superior performance over other frequently used methods. 

## Citation

Please cite the GAGE paper when using this open-source  package. This will help the project and our team:

Luo W, Friedman M, etc. GAGE: generally applicable gene set enrichment for pathway analysis. BMC Bioinformatics, 2009, 10, pp. 161, <a href=https://doi.org/10.1186/1471-2105-10-161>doi: 10.1186/1471-2105-10-161</a>

## Installation

``` r
# install from BioConductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("gage")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("datapplab/gage")
```

## Quick start with demo data

``` r
#preparation
library(gage)
data(gse16873)
hn=(1:6)*2-1
dcis=(1:6)*2

#KEGG pathway analysis
data(kegg.gs)
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis)

#GO term analysis, separate BP, MF and CC categories, need the gageData pacakge
library(gageData)
data(go.sets.hs)
data(go.subs.hs)
gse16873.bp.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$BP], ref = hn, samp = dcis)
gse16873.mf.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$MF], ref = hn, samp = dcis)
gse16873.cc.p <- gage(gse16873, gsets = go.sets.hs[go.subs.hs$CC], ref = hn, samp = dcis)
```

## More information

Please check the <a href=https://bioconductor.org/packages/gage/>BioC page</a> for tutorials and extra documentations. Thank you for your interest.

