#' @title regionSet
#'
#' @description Stores .cov and .cnv file (BED inspired file) into a RaggedExperiment object. Will use genomic locations from the .cnv file to find the start and end points of the CNVs from the .cov file and store these as two lists. The assay of the RE object will store the cov file and a list of nested dataframes called "CNVlist" within the metadata will store all the cnv dataframes.
#'
#' @param cnv List of CNVs in a dataframe containing CNVs from detection algorithms/ pipelines. Additional columns such as BATCH and GENE can also be used for better plotting.
#' @param cov cov file from WES platform/ sequencer kit or if WGS regular intervals. Stored
#' as a dataframe - genomic locations as rows and samples as columns. It is recommended to normalize Read Depth by library size for fairer comparisons.
#' @param col colData - A string which is used in RaggedExperiment creation. Default is "experiment".
#'
#' @return RaggedExperiment object with two lists: Start sites and End sites.The cov file will be stored as the main assay of the RE object, and a new list will appear in the metadata to contain all cnv dataframes.
#'
#' @import RaggedExperiment
#' @import GenomicRanges
#' @export
#'
#' @examples
#' data("test_cnv")
#' data("test_cov")
#' SARC <- regionSet(cnv = test_cnv, cov = test_cov)
regionSet <- function(cnv, cov, col = "experiment"){

  if (missing(cnv)) stop('cnv is missing. Add cnv dataframe. Ideally the most recently created cnv file should be used.')

  if (missing(cov)) stop('cov is missing. Add a cov dataframe. This is a coverage file consisting of genomic locations and read depths (ideally normalised by library size) for each sample in the WES/WGS cohort.')

  #convert cov to a grange object and then a RaggedExperiment

  gr <- makeGRangesFromDataFrame(cov, keep.extra.columns = TRUE)

  gr <- sortSeqlevels(gr)

  gr <- sort(gr)

  ragexp <- RaggedExperiment(colData = col, gr)

  #create objects for storage

  cov.mini <- list()

  start.mini <- vector()

  end.mini <- vector()

  for (i in seq_len(nrow(cnv))) {

    #Break cov into chromosomes - for speed

    cov.mini[[i]] <- ragexp@assays@unlistData[seqnames(ragexp@assays@unlistData) == cnv$CHROM[i]]

    #Make iteration objects

    c <- cov.mini[[i]]

    s <- cnv$START[i]

    e <- cnv$END[i]

    #Get ith START from cov
    s.min <- which.min(abs(start(c@ranges)-s))

    start.mini[i] <- s.min

    #Get ith END from cov
    e.min <- which.min(abs(end(c@ranges)-e))

    end.mini[i] <- e.min
  }

  #create list to keep cnv dataframes
  cnvlist <- list()

  cnvlist[["Original"]] <- cnv

  metadata(ragexp)[["CNVlist"]] <- cnvlist

  #store start and end regions in RE

  metadata(ragexp)[["REGIONSTART"]] <- as.numeric(start.mini)

  metadata(ragexp)[["REGIONEND"]] <- as.numeric(end.mini)

  return(ragexp)

}
