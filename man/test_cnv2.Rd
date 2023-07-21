\name{test_cnv2}
\alias{test_cnv2}
\docType{data}
\title{
test cnv 2
}
\description{
This is not a traditional bed file, but shares the purpose of a bed file during CNV analysis - which is to store a list of CNVs detected from CNV detection algorithms/ pipelines. This files shares the first three columns of a traditional bed file (chrom, chromstart, chromend), and in addition has several more columns e.g. sample, type,value/ score and tool. The column names must match: SAMPLE, CHROM, START, END, TYPE, VALUE - but only the order of the first three columns must remain the same. In addition, columns for the batch and gene are included to improve visualisation.
}
\usage{data("test_cnv2")}
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
  }
}
\details{
Several annonymised whole exome sequencing samples which have been put through a Copy Number Variaiton pipleine and had their CNVs recorded in a way which resembles a bed file. This file has additional columns too.
}
\source{
This is annonmysed patient data from the Newcastle Mitochondrial Research group, from the Wellcome Center for Mitochondrial research, Newcastle Universiry, UK.
}
\references{
Taylor, Robert W., et al. "Use of whole-exome sequencing to determine the genetic basis of multiple mitochondrial respiratory chain complex deficiencies." Jama 312.1 (2014): 68-77.
}
\value{
Dataframe listing CNVs from a WES CNV detection pipeline. Contains extra columns for plot specificity.
}
\examples{
data(test_cnv2)
## maybe str(test_cnv2) ; plot(test_cnv2) ...
}
\keyword{datasets}
