#' @title splitHelp
#'
#' @description Internal function to aid regionSplit to split the large cov file into multiple smaller cov files - one for each CNV in the cnv file.
#'
#' @param cov cov file from WES platform/ sequencer kit or if WGS regular intervals. Stored
#' as a dataframe - genomic locations as rows and samples as columns.
#'
#' @param x Genomic start site.
#'
#' @param x Genomic end site.
#'
#' @return Helps regionSplit to return smaller cov dataframes, one for each CNV.
#'
#' @noRd

splitHelp <- function(cov,x,y){

  if (missing(cov)) stop('cov is missing. Add a cov dataframe. This is a coverage file consisting of genomic locations and read depths for each sample in the WES/WGS cohort.')

  if (missing(x)) stop('x is missing. Add an integer representing the start position of a CNV.')

  if (missing(y)) stop('y is missing. Add an integer representing the end position of a CNV.')

  #function for mapply in SARC::regionSplit

  z <- cov[x:y,]

  #return smaller cov file/ dataframe

  return(z)

}
