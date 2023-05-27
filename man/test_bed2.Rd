\name{test_bed2}
\alias{test_bed2}
\docType{data}
\title{
test bed 2.
}
\description{
An example data file for testing further functions allowed by the SARC R package.
Use of BATCH, GENE and VCF columns can lead to more specific and better images
for viewing and publications. This file gives examples as to how users should organise their data files for processing by SARC.
}
\usage{data("test_bed2")}
\format{
  A data frame with 15 observations on the following 10 variables.
  \describe{
    \item{\code{SAMPLE}}{a character vector}
    \item{\code{CHROM}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
    \item{\code{START}}{a numeric vector}
    \item{\code{END}}{a numeric vector}
    \item{\code{TYPE}}{a factor with levels \code{DEL} \code{DUP}}
    \item{\code{VALUE}}{a numeric vector}
    \item{\code{TOOL}}{a factor with levels \code{clinCNV} \code{cn.mops} \code{exomeDepth}}
    \item{\code{BATCH}}{a character vector}
    \item{\code{GENE}}{a character vector}
    \item{\code{VCF}}{a character vector}
  }
}
\details{
Our pipeline used clinCNV, cn.mops and exomedepth. Any CNV detection tool can be used.
We only kept entries from chromosome 1 to keep the test data small.
}
\value{A dataframe which represents a BED file which is the most common output type in CNV detection pipelines/ algorithms. This object will list the key information for each CNV, including the sample, chromosome, start and end sites, type of cnv and the cnv value. Further columns (BATCH, GENE, VCF) were added
here to show extra functions.}

\examples{
data(test_bed2)
## maybe str(test_bed2) ; plot(test_bed2) ...
}
\keyword{datasets}
