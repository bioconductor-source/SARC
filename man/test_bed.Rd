\name{test_bed}
\alias{test_bed}
\docType{data}
\title{
test .bed file.
}
\description{
This is a generic .bed file generally used as a common representation of CNVs detected from CNV detection pipelines. Columns should be: Sample IDs, Chromosome, Start, End, Type of CNV, CNV value (must be represented as an interger e.g. 0=HOM DELETION, 1=HET DELETION, 2=Normal, 3=HET DUPLICATION, 4/>4=HOM DUPLICATION), and TOOL. Additional columns such as gene of interest can also be added.
}
\usage{data("test_bed")}
\format{
  A data frame with 15 observations on the following 7 variables.
  \describe{
    \item{\code{SAMPLE}}{a character vector}
    \item{\code{CHROM}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrX} \code{chrY}}
    \item{\code{START}}{a numeric vector}
    \item{\code{END}}{a numeric vector}
    \item{\code{TYPE}}{a factor with levels \code{DEL} \code{DUP}}
    \item{\code{VALUE}}{a numeric vector}
    \item{\code{TOOL}}{a factor with levels \code{clinCNV} \code{cn.mops} \code{exomeDepth}}
  }
}
\details{
Our pipeline used clinCNV, cn.mops and exomedepth. Any CNV detection tool can be used.
We only kept entries from chromosome 1 to keep the test data small.
}

\value{A dataframe which represents a BED file which is the most common output type in CNV detection pipelines/ algorithms. This object will list the key information for each CNV, including the sample, chromosome, start and end sites, type of cnv and the cnv value.}

\examples{
data(test_bed)
## maybe str(test_bed) ; plot(test_bed) ...
}
\keyword{datasets}
