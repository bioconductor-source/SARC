#' altPos
#'
#' @description Internal function for plotCovPrep. Here the CNV region is padded by x rows to make it easier to distinguish a large CNV variations in the detected regions.
#'
#' @param start Start site of a CNV.
#' @param end End site of a CNV.
#' @param n1 How many rows to pad the start site of each CNV - dafault is 0.
#' @param n2 How many rows to pad the end site of each CNV - dafault is 0.
#'
#' @return Smaller cov files with padded start and end reads.
#' @noRd
altPos <- function(start, end, n1=0, n2=0){

  if (missing(start)) stop('start is missing')

  if (missing(end)) stop('end is missing.')

  #start padding
  start <- start - n1

  #end padding
  end <- end + n2

  #create new start and end positions for each CNV
  pos <- cbind(start, end)

  #return new positions
  return(pos)

}
