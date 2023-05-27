#' prepBase
#'
#' @description Internal function for prepAnova.This function avoids errors in the start - end site differences. It slightly increases the size of the CNVs if the CNV is very small - and only covers < 3 reads / rows of the cov file.
#'
#' @param start Start site for CNV.
#' @param end End site for CNV.
#'
#' @noRd
prepBase <- function(start, end){

  if (missing(start)) stop('start is missing')

  if (missing(end)) stop('end is missing.')

  #if CNV region is very short, make it slightly longer

  if (isTRUE(end - start == 0)) {

    start <- start - 1

    end <- end + 1

  } else if (isTRUE(end - start == 1)) {

    start <- start - 1

  }

  prepped <- cbind(start, end)

  return(prepped)

}
