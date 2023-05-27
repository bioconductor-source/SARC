\name{test_cov}
\alias{test_cov}
\docType{data}
\title{
test .cov file.
}
\description{
The cov/ coverage file can be generated using bcftools, samtools among other genomic tools. It uses the WES platform .bed file from the chemistry of the Illumina DNA Rapid test. The first column should be ID (same as rownames), Chromosome, Start, End and the rest of the columns represent the samples. Intergers will represent the number of reads from each samples .bam file.
}
\usage{data("test_cov")}
\format{
  A data frame with 34987 observations on the following 24 variables.
  \describe{
    \item{\code{ID}}{a numeric vector}
    \item{\code{CHROM}}{a factor with levels \code{chr1} \code{chr10} \code{chr11} \code{chr12} \code{chr13} \code{chr14} \code{chr15} \code{chr16} \code{chr17} \code{chr18} \code{chr19} \code{chr2} \code{chr20} \code{chr21} \code{chr22} \code{chr3} \code{chr4} \code{chr5} \code{chr6} \code{chr7} \code{chr8} \code{chr9} \code{chrM} \code{chrX} \code{chrY}}
    \item{\code{START}}{a numeric vector}
    \item{\code{END}}{a numeric vector}
    \item{\code{Sample_616}}{a numeric vector}
    \item{\code{Sample_617}}{a numeric vector}
    \item{\code{Sample_618}}{a numeric vector}
    \item{\code{Sample_619}}{a numeric vector}
    \item{\code{Sample_620}}{a numeric vector}
    \item{\code{Sample_621}}{a numeric vector}
    \item{\code{Sample_622}}{a numeric vector}
    \item{\code{Sample_623}}{a numeric vector}
    \item{\code{Sample_624}}{a numeric vector}
    \item{\code{Sample_625}}{a numeric vector}
    \item{\code{Sample_626}}{a numeric vector}
    \item{\code{Sample_627}}{a numeric vector}
    \item{\code{Sample_628}}{a numeric vector}
    \item{\code{Sample_629}}{a numeric vector}
    \item{\code{Sample_630}}{a numeric vector}
    \item{\code{Sample_635}}{a numeric vector}
    \item{\code{Sample_636}}{a numeric vector}
    \item{\code{Sample_637}}{a numeric vector}
    \item{\code{Sample_638}}{a numeric vector}
    \item{\code{Sample_639}}{a numeric vector}
  }
}
\details{
Only chromosome 1 was kept to keep the test data small.
Samples are annonymised patient samples used for CNV detection.
}

\value{
A dataframe which represents a COV file. This is one of the most common input file types for CNV detection algorithms. This is generated from BAM files from WES/ WGS experiments. Regions of the genome have their read-depths/ coverage counted for each sample.
}

\examples{
data(test_cov)
## maybe str(test_cov) ; plot(test_cov) ...
}
\keyword{datasets}
