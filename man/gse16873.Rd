\name{gse16873}
\Rdversion{1.1}
\alias{gse16873}
\docType{data}
\title{
GSE16873: a breast cancer microarray dataset
}
\description{
GSE16873 is a breast cancer study (Emery et al, 2009) downloaded from
Gene Expression Omnibus (GEO). GSE16873 covers twelve patient cases,
each with HN (histologically
 normal), ADH (ductal hyperplasia), and DCIS (ductal carcinoma in situ)
 RMA samples. Due to the size limit of this package, we split this
 GSE16873 into two halves, each including 6 patients with their HN and
 DCIS but not ADH tissue types. This gage package only
includes the first half dataset for 6 patients as this example dataset
gse16873. Most of our demo analyses are done on the first half dataset,
except for the advanced analysis where we use both halves datasets with
all 12 patients. Details section below gives more information.
}

\usage{data(gse16873)}
\details{
Raw data for these two half datasets were
 processed separately using two different methods, FARMS
and RMA, respectively to generate the non-biological data
heterogeneity. The first half dataset is named as gse16873, the second
half dataset named gse16873.2. We also have the full dataset,
gse16873.full, which includes all HN, ADH and DCIS samples of all 12
patients, processed together using FARMS. The second half dataset
plus the full dataset and the original BMP6 dataset
used in GAGE paper and earlier versions of gage is supplied
with another package, gageData.
}
\source{
  GEO Dataset GSE16873:
  <URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16873>
}
\references{
Emery LA, Tripathi A, King C, Kavanah M, Mendez J, Stone MD, de las
Morenas A, Sebastiani P, Rosenberg CL: Early dysregulation of cell
adhesion and extracellular matrix pathways in breast cancer
progression. Am J Pathol 2009, 175:1292-302.
}
\examples{
data(gse16873)
#check the heterogenity of the two half datasets
boxplot(data.frame(gse16873))

#column/smaple names
cn=colnames(gse16873)
hn=grep('HN',cn, ignore.case =TRUE)
adh=grep('ADH',cn, ignore.case =TRUE)
dcis=grep('DCIS',cn, ignore.case =TRUE)
print(hn)
print(dcis)

data(kegg.gs)
lapply(kegg.gs[1:3],head)
head(rownames(gse16873))
gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs,
    ref = hn, samp = dcis)
}
\keyword{datasets}
