#' setupPlotStep2
#'
#' @param covprepped covprepped List of dataframes, one for each CNV. This list should be in metadata after bring created by SARC::plotCovPrep.
#' @param sps1 List of dataframes which include genes and exons for each CNV region. Output from SARC::setupPlotStep1.
#' @param START Column which is required to setkey from data.table.
#'
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @noRd
setupPlotStep2 <- function(covprepped, sps1, START = "START"){

  if (missing(covprepped)) stop('covprepped is missing. Add a list of small cov dataframes that have been prepared for plotting. Should be created from SARC::plotCovPrep.')

  if (missing(sps1)) stop('sps1 is missing. This is output from setupPlotStep1.')

  #set the START column as a key

  pl <- lapply(covprepped, function(x){

    x <- as.data.frame(x)

    names(x)[1:3] <- c("CHROM", "START", "END")

    y <- data.table(x)

    setkey(y, START)

    return(y)

  })

  #alter names of the dataframes

  yl <- lapply(sps1, function(x){

    names(x)[1:3] <- c("CHROM", "START", "END")

    y <- data.table(x)

    setkey(y, START)

    return(y)

  })

  #re-organise the dataframes by nucleotide start positions

  xl <- list()

  for (i in seq_len(length(yl))) {

    xl[[i]] <- yl[[i]][pl[[i]], roll='nearest', on= "START"]

    xl[[i]]$LOC <- paste0(xl[[i]]$START, "_", xl[[i]]$GENE_EXON)

    xl[[i]] <- xl[[i]][,-c(1:14)]

  }

  #apply setupPlotStep2fun

  ml <- lapply(xl, function(x) setupPlotStep2fun(x))

  #return new list

  return(ml)
}
