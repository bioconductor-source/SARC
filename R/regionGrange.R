#' regionGrange
#'
#' @description Internal function used for regionGrangeMake. This function will create a grange object for each region where a CNV was detected.
#'
#' @param minicov List of dataframes, each dataframe is small cov file which ranges the detected CNV. This will be a list stored as metadata after SARC::plotCovPrep.
#' @param range Number of rows (genomic ranges) to group together as one entry for grange. As some WES platforms will have very short genomic start-end coordinates, this is useful to extract the exon/ gene data. Default is 10.
#' @param gap The gap between the start of one grange entry to the next entry. Default is 10.
#'
#' @noRd
#'
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
regionGrange <- function(minicov, range=10, gap=10){

  if (missing(minicov)) stop('minicov is missing. Add list of small cov file, postentially with padding. Should be after SARC::plotCovPrep.')

  #make sure the CNV/ reigon is only from a single chromosome

  x <- minicov

  x <- as.data.frame(x)

  chrom <- as.character(x$seqnames[1])

  x <- x[which(x$seqnames == chrom),]

  #create grange parameters

  seqGrange <- seq(1, nrow(x), range)

  seqGrange <- as.data.frame(cbind(seqGrange, seqGrange+gap))

  colnames(seqGrange) <- c("start", "end")

  seqGrange[,2][nrow(seqGrange)] <- nrow(x)

  #loop add start and end sites to the dataframe of start and end sites

  for (i in seq_len(nrow(seqGrange))) {

    seqGrange$start[i] <- x$start[seqGrange$start[i]]

    seqGrange$end[i] <- x$end[seqGrange$end[i]]

  }

  #create grange object for each row of the dataframe

  GrangeList <- list()

  for (i in seq_len(nrow(seqGrange))) {

    GrangeList[[i]] <- GRanges(seqnames = chrom,

                               strand = "+",

                               ranges = IRanges(start = seqGrange[[1]][i],

                                                end = seqGrange[[2]][i]))

    names(GrangeList)[i] <- paste0(chrom, ":", seqGrange[[1]][i], ":",

                                   seqGrange[[2]][i])
  }

  #return a list of grange objects

  return(GrangeList)

}
