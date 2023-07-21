\name{test_cov}
\alias{test_cov}
\docType{data}
\title{
test coverage file
}
\description{
This file is a large matrix of raw read depth from whole exome sequencing BAM files. The first four columns are: ID, chromosome, start position and end position.Following on, each sample has its own column which gives an integer for each position of DNA found in the whole exome experiment. Samples are annoymised. Read depth has been normalised by total library sizes via RPM (reads per million) calculations.
}
\usage{data("test_cov")}
\format{
  A data frame with 34987 observations on the following 24 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{CHROM}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrM} \code{chrX} \code{chrY}}
    \item{\code{START}}{a numeric vector}
    \item{\code{END}}{a numeric vector}
    \item{\code{SampleA}}{a numeric vector}
    \item{\code{SampleB}}{a numeric vector}
    \item{\code{SampleC}}{a numeric vector}
    \item{\code{SampleD}}{a numeric vector}
    \item{\code{SampleE}}{a numeric vector}
    \item{\code{SampleF}}{a numeric vector}
    \item{\code{SampleG}}{a numeric vector}
    \item{\code{SampleH}}{a numeric vector}
    \item{\code{SampleI}}{a numeric vector}
    \item{\code{SampleJ}}{a numeric vector}
    \item{\code{SampleK}}{a numeric vector}
    \item{\code{SampleL}}{a numeric vector}
    \item{\code{SampleM}}{a numeric vector}
    \item{\code{SampleN}}{a numeric vector}
    \item{\code{SampleO}}{a numeric vector}
    \item{\code{SampleP}}{a numeric vector}
    \item{\code{SampleQ}}{a numeric vector}
    \item{\code{SampleR}}{a numeric vector}
    \item{\code{SampleS}}{a numeric vector}
    \item{\code{SampleT}}{a numeric vector}
  }
}
\details{
Several annonymised whole exome sequencing samples had been converted to a bam file and had a read counting tool applied.Only values from chromosome 1 are present to save space and time in examples.
}
\source{
This is annonmysed patient data from the Newcastle Mitochondrial Research group, from the Wellcome Center for Mitochondrial research, Newcastle Universiry, UK.
}
\references{
Taylor, Robert W., et al. "Use of whole-exome sequencing to determine the genetic basis of multiple mitochondrial respiratory chain complex deficiencies." Jama 312.1 (2014): 68-77.
}
\value{
Dataframe of normalised coverage from WES samples.
}
\examples{
data(test_cov)
## maybe str(test_cov) ; plot(test_cov) ...
}
\keyword{datasets}
